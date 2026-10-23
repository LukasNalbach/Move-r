/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_TIMEOUT_RUNNER_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_TIMEOUT_RUNNER_HPP

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <ctime>
#include <optional>
#include <string>

#include <poll.h>
#include <sys/wait.h>
#include <unistd.h>

#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/memory.hpp>

// the result of a query: its number of occurrences and its time in ns
struct timed_result_t {
    uint64_t count;
    uint64_t elapsed_ns;
};

/**
 * @brief reads exactly n bytes from a file descriptor
 * @return false on EOF (the other end of the pipe is closed) or on an error
 */
inline bool read_full(int32_t fd, void* buf, uint64_t n)
{
    char* p = static_cast<char*>(buf);

    while (n > 0) {
        int64_t r = read(fd, p, n);
        if (r > 0) { p += r; n -= (uint64_t) r; }
        else if (r == 0) return false;
        else if (errno == EINTR) continue;
        else return false;
    }

    return true;
}

/**
 * @brief writes exactly n bytes to a file descriptor
 * @return false if the other end of the pipe is closed or on an error
 */
inline bool write_full(int32_t fd, const void* buf, uint64_t n)
{
    const char* p = static_cast<const char*>(buf);

    while (n > 0) {
        int64_t w = write(fd, p, n);
        if (w > 0) { p += w; n -= (uint64_t) w; }
        else if (w < 0 && errno == EINTR) continue;
        else return false;
    }

    return true;
}

/**
 * @brief waits up to timeout_ns ns until a file descriptor can be read from
 * @return false on a timeout or an error
 */
inline bool wait_readable(int32_t fd, uint64_t timeout_ns)
{
    struct pollfd pfd { fd, POLLIN, 0 };
    constexpr uint64_t ns_per_second = 1'000'000'000;
    struct timespec ts { (time_t) (timeout_ns / ns_per_second), (int64_t) (timeout_ns % ns_per_second) };
    int32_t r;
    do { r = ppoll(&pfd, 1, &ts, nullptr); } while (r < 0 && errno == EINTR);
    return r > 0 && (pfd.revents & POLLIN);
}

/**
 * @brief runs queries in a child process (created with fork(), so it shares the loaded index) and kills it if a
 *        query takes longer than a timeout. A query is a request of request_size bytes, sent through a pipe; the child
 *        answers with the timed_result_t of query(request).
 * @tparam query_fnc_t called in the child as query(request) with a std::string of request_size bytes; returns the
 *         timed_result_t of the query (the time measured in the child)
 */
template <typename query_fnc_t>
class timeout_runner {
    query_fnc_t query;
    uint64_t request_size; // number of bytes of a request
    pid_t child = -1; // process id of the child
    int32_t to_child = -1; // write end of the pipe of the requests
    int32_t from_child = -1; // read end of the pipe of the results

    // creates the child and the pipes
    void spawn()
    {
        int32_t down[2], up[2]; // down: requests, up: results
        if (pipe(down) != 0 || pipe(up) != 0) { perror("pipe"); std::exit(1); }
        pid_t pid = fork();
        if (pid < 0) { perror("fork"); std::exit(1); }

        if (pid == 0) {
            // the child answers requests until the parent closes the pipe or kills it
            close(down[1]); close(up[0]);
            std::string request;
            resize_no_init(request, request_size);

            while (read_full(down[0], request.data(), request_size)) {
                timed_result_t r = query(request);
                if (!write_full(up[1], &r, sizeof(r))) break;
            }

            _exit(0);
        }

        close(down[0]); close(up[1]);
        child = pid; to_child = down[1]; from_child = up[0];
    }

    // kills the child and waits for it
    void reap()
    {
        if (child <= 0) return;
        kill(child, SIGKILL);
        int32_t st;
        while (waitpid(child, &st, 0) < 0 && errno == EINTR) {}
        close(to_child); close(from_child);
        child = -1; to_child = from_child = -1;
    }

public:
    timeout_runner(query_fnc_t query, uint64_t request_size) : query(std::move(query)), request_size(request_size)
    {
        spawn();
    }

    ~timeout_runner() { reap(); }

    /**
     * @brief runs a query in the child and waits up to timeout_ns ns for its result
     * @param request the request (request_size bytes)
     * @param timeout_ns the timeout in ns
     * @return the result, or std::nullopt if the query timed out (then the child is killed and a new one is created)
     */
    std::optional<timed_result_t> run(const std::string& request, uint64_t timeout_ns)
    {
        if (write_full(to_child, request.data(), request_size) && wait_readable(from_child, timeout_ns)) {
            timed_result_t r;
            if (read_full(from_child, &r, sizeof(r))) return r;
        }

        reap();
        spawn();
        return std::nullopt;
    }
};

/**
 * @brief writes a pattern file (a pizza&chili header and the patterns) and removes it from the page cache
 * @param path path of the pattern file
 * @param patterns the concatenated patterns
 * @param num_patterns number of patterns
 * @param m pattern length
 * @param text_name name of the text (in the header)
 * @param forbidden the symbols the patterns do not contain (in the header)
 */
inline void write_pattern_file(const std::string& path, const std::string& patterns, uint64_t num_patterns, uint64_t m,
    const std::string& text_name, const std::string& forbidden)
{
    direct_ofstream file(path, std::ios::binary);
    file << "# number=" << num_patterns << " length=" << m << " file=" << text_name << " forbidden=" << forbidden << "\n";
    file.write(patterns.data(), patterns.size());
    file.close();
    drop_from_page_cache(path, true);
}

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_TIMEOUT_RUNNER_HPP
