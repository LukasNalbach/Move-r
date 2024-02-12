#pragma once

#include <array>
#include <vector>
#include <iostream>
#include <cstring>

/**
 * @brief variable-width interleaved vectors (widths are fixed to whole bytes)
 * @tparam uint_t unsigned integer type
 */
template <typename uint_t = uint32_t, uint8_t num_vectors = 8>
class interleaved_vectors {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);
    static_assert(num_vectors > 0);

    protected:
    uint64_t size_vectors = 0; // size of each stored vector
    uint64_t capacity_vectors = 0; // capacity of each stored vector
    uint64_t width_entry = 0; // sum of the widths of all vectors

    // [0..(capacity_vectors+1)*num_vectors-1] vector storing the interleaved vectors
    std::vector<char> data_vectors;

    // [0..num_vectors-1] or shorter; widths of the stored vectors; widths[i] = width of vector i
    std::array<uint64_t,num_vectors> widths;

    /** @brief [0..num_vectors-1] pointers to the first entries of each vector; bases[vec] = base
     *        of the first entry of the vector with index vec */
    std::array<char*,num_vectors> bases;

    /** @brief [0..num_vectors-1] masks that are used to mask off data of other vector entries when
     *        accessing a vector */
    std::array<uint_t,num_vectors> masks;

    /**
     * @brief initializes the interleaved_vectors with the vector-widths stored in widths
     * @param widths vector containing the widths (in bytes) of the interleaved arrays
     */
    void initialize(std::array<uint8_t,num_vectors> widths = {sizeof(uint_t)}) {
        size_vectors = 0;
        capacity_vectors = 0;
        width_entry = 0;
        
        data_vectors.clear();
        data_vectors.shrink_to_fit();

        for (uint8_t i=0; i<num_vectors; i++) {
            if (widths[i] == 0) {
                masks[i] = 0;
                this->widths[i] = 0;
            } else {
                this->widths[i] = widths[i];
                width_entry += widths[i];
                masks[i] = std::numeric_limits<uint_t>::max()>>(8*(sizeof(uint_t)-widths[i]));
            }
        }

        reserve(2);
    }
    
    /**
     * @brief builds the bases array
     * @param data pointer to the data storing the interleaved vectors
     */
    void set_bases(char* data) {
        bases[0] = data;

        for (uint8_t i=1; i<num_vectors; i++) {
            if (widths[i] == 0) {
                bases[i] = NULL;
            } else {
                bases[i] = bases[i-1] + widths[i-1];
            }
        }
    }

    /**
     * @brief copies another interleaved_vectors object into this object
     * @param other another interleaved_vectors object
     */
    void copy_from_other(const interleaved_vectors& other) {
        size_vectors = other.size_vectors;
        capacity_vectors = other.capacity_vectors;
        width_entry = other.width_entry;
        data_vectors = other.data_vectors;
        widths = other.widths;
        masks = other.masks;

        set_bases(&data_vectors[0]);
    }

    /**
     * @brief moves another interleaved_vectors object into this object
     * @param other another interleaved_vectors object
     */
    void move_from_other(interleaved_vectors&& other) {
        size_vectors = other.size_vectors;
        capacity_vectors = other.capacity_vectors;
        width_entry = other.width_entry;

        data_vectors = std::move(other.data_vectors);
        widths = other.widths;
        bases = other.bases;
        masks = other.masks;

        other.initialize();
    }
    
    public:
    interleaved_vectors() {initialize();}
    interleaved_vectors(interleaved_vectors&& other) {move_from_other(std::move(other));}
    interleaved_vectors(const interleaved_vectors& other) {copy_from_other(other);}
    interleaved_vectors& operator=(interleaved_vectors&& other) {move_from_other(std::move(other));return *this;}
    interleaved_vectors& operator=(const interleaved_vectors& other) {copy_from_other(other);return *this;}

    ~interleaved_vectors() {
        for (uint8_t i=0; i<num_vectors; i++) {
            bases[i] = NULL;
        }
    }

    /**
     * @brief Construct a new interleaved_vectors 
     * @param widths vector containing the widths (in bytes) of the interleaved arrays
     */
    interleaved_vectors(std::array<uint8_t,num_vectors> widths) {
        initialize(widths);
    }

    /**
     * @brief returns the size of each stored vector
     * @return the size of each stored vector
     */
    uint64_t size() const {
        return size_vectors;
    }

    /**
     * @brief returns whether the interleaved vectors are empty
     * @return whether the interleaved vectors are empty
     */
    bool empty() const {
        return size_vectors == 0;
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const {
        return
            2*sizeof(uint64_t)+1+ // variables
            num_vectors*sizeof(uint64_t)+ // widths
            2*num_vectors*sizeof(uint_t)+ // masks
            size_vectors*width_entry; // data_vectors
    }

    /**
     * @brief returns total width (number of bytes) per entry, that is the (sum of all widths)
     * @return number of bytes per entry
     */
    uint8_t bytes_per_entry() const {
        return this->width_entry;
    }

    /**
     * @brief returns the width in bytes of the vector with index vec
     * @param vec vector index
     * @return its witdth in bytes
     */
    uint8_t width(uint8_t vec) const {
        return widths[vec];
    }

    /**
     * @brief returns a pointer to the data of the interleved vectors
     * @return pointer to the data of the interleved vectors
     */
    char* data() const {
        return bases[0];
    }

    /**
     * @brief returns the i-th entry of the vector with index 0
     * @param i entry index (0 <= i < size_vectors)
     * @return i-th entry of the vector with index 0 
     */
    uint_t operator[](uint_t i) const {
        return get<0>(i);
    }

    /**
     * @brief reserves capacity entries in all stored vectors; if capacity is smaller than the
     *        current capacity, nothing happens (use shrink_to_fit() to lower the capacity)
     * @param capacity capacity
     * @param num_threads number of threads to use when copying entries (default: 1)
     */
    void reserve(uint64_t capacity, uint8_t num_threads = 1) {
        if (capacity_vectors < capacity) {
            std::vector<char> new_data_vectors;
            no_init_resize(new_data_vectors,(capacity+1)*width_entry);

            #pragma omp parallel for num_threads(num_threads)
            for (uint64_t i=0; i<size_vectors*width_entry; i++) {
                new_data_vectors[i] = data_vectors[i];
            }

            std::memset(&new_data_vectors[size_vectors*width_entry],0,width_entry);
            std::memset(&new_data_vectors[capacity*width_entry],0,width_entry);
            std::swap(data_vectors,new_data_vectors);
            set_bases(&data_vectors[0]);
            capacity_vectors = capacity;
        }
    }

    /**
     * @brief resizes all stored vectors to size and initializes new entries to 0
     * @param size size
     * @param num_threads number of threads to use when copying entries and initializing new entries
     *                    to 0 (default: 1)
     */
    void resize(uint64_t size, uint8_t num_threads = 1) {
        uint64_t size_vectors_old = size_vectors;

        if (capacity_vectors < size) {
            reserve(size,num_threads);
        }
        
        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i=size_vectors_old*width_entry; i<size*width_entry; i++) {
            data_vectors[i] = uchar_to_char(0);
        }

        size_vectors = size;
    }

    /**
     * @brief resizes all stored vectors to size without initializing new entries to 0
     * @param size size
     * @param num_threads number of threads to use when copying entries (default: 1)
     */
    void resize_no_init(uint64_t size, uint8_t num_threads = 1) {
        if (capacity_vectors < size) {
            reserve(size,num_threads);
        }

        size_vectors = size;
    }

    /**
     * @brief resizes all stored vectors to size 0
     */
    void clear() {
        resize(0);
    }

    /**
     * @brief shrinks all stored vectors to their size
     */
    void shrink_to_fit() {
        if (size_vectors < capacity_vectors) {
            capacity_vectors = std::max<uint64_t>(2,size_vectors);
            data_vectors.resize((capacity_vectors+1)*width_entry);
            data_vectors.shrink_to_fit();
            set_bases(&data_vectors[0]);
        }
    }

    /**
     * @brief appends a tuple of vec values to the end of the interleaved vectors
     * @param vals tuple of vec values
     */
    template <typename... Ts>
    void emplace_back(std::tuple<Ts...>&& vals) {
        static_assert(std::tuple_size<std::tuple<Ts...>>::value <= num_vectors);
        
        if (size_vectors == capacity_vectors) {
            reserve(1.5*size_vectors);
        }

        for_constexpr<0,sizeof...(Ts),1>([this,&vals](auto vec){
            set<vec>(size_vectors,vals[vec]);
        });

        size_vectors++;
    }

    /**
     * @brief appends a tuple of vec values to the end of the interleaved vectors
     * @param values tuple of vec values
     */
    template <typename... Ts>
    void push_back(std::tuple<Ts...> values) {
        emplace_back<Ts...>(std::move(values));
    }

    /**
     * @brief unsafely appends a tuple of vec values to the end of the interleaved vectors (if
     *        there exists a vec \in [0,sizeof...(Ts...)) s.t. sizeof(Ts...[vec]) > widths[vec],
     *        this is unsafe, because subsequent entries are overwritten)
     * @param vals tuple of vec values
     */
    template <typename... Ts>
    void emplace_back_unsafe(std::tuple<Ts...>&& vals) {
        static_assert(std::tuple_size<std::tuple<Ts...>>::value <= num_vectors);

        if (size_vectors == capacity_vectors) {
            reserve(1.5*size_vectors);
        }

        for_constexpr<0,sizeof...(Ts),1>([this,&vals](auto vec){
            set_unsafe<vec,std::tuple_element_t<vec,std::tuple<Ts...>>>(size_vectors,std::get<vec>(vals));
        });

        size_vectors++;
    }

    /**
     * @brief appends a tuple of vec values to the end of the interleaved vectors
     * @param vals tuple of vec values
     */
    template <typename... Ts>
    void push_back_unsafe(std::tuple<Ts...> vals) {
        emplace_back_unsafe<Ts...>(std::move(vals));
    }

    /**
     * @brief sets the i-th entry in the vector with index vec to v
     * @tparam vec vector index (0 <= vec < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec>
    inline void set(uint_t i, uint_t v) {
        static_assert(vec < num_vectors);

        for (uint64_t byte=0; byte<widths[vec]; byte++) {
            *reinterpret_cast<char*>(bases[vec] + (i * width_entry) + byte) =
            *(reinterpret_cast<char*>(&v) + byte);
        }
    }

    /**
     * @brief plainly stores the value v (of type T) at the memory location of the i-th entry in the vector
     * with index vec (if sizeof(T) > widths[vec], this is unsafe, since subsequnet entries are overwritten)
     * @tparam vec vector index (0 <= vec < num_vectors)
     * @tparam T type to treat the i-th entry in the vector with index vec as
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec, typename T>
    inline void set_unsafe(uint_t i, T v) {
        static_assert(vec < num_vectors);
        *reinterpret_cast<T*>(bases[vec] + i * width_entry) = v;
    }

    /**
     * @brief returns the i-th entry in the vector with index vec
     * @tparam vec vector index (0 <= vec < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @return value
     */
    template <uint8_t vec>
    inline uint_t get(uint_t i) const {
        static_assert(vec < num_vectors);
        return *reinterpret_cast<uint_t*>(bases[vec] + i * width_entry) & masks[vec];
    }

    /**
     * @brief interprets the memory location of the i-th entry in the vector with index vec as an object of type
     * T and returns it (if sizeof(T) > widths[vec], the returned value contains data from subsequent entries)
     * @tparam vec vector index (0 <= vec < num_vectors)
     * @tparam uint_t type to treat the i-th entry in the vector with index vec as
     * @param i entry index (0 <= i < size_vectors)
     * @return value
     */
    template <uint8_t vec, typename T>
    inline uint_t get_unsafe(uint_t i) const {
        static_assert(vec < num_vectors);
        return *reinterpret_cast<T*>(bases[vec] + i * width_entry);
    }

    /**
     * @brief reinterpret the memory at data as interleaved vectors of size size; do not perform any operations that
     * may change the size or the capacity of the interleaved vectors after using this method
     * @param data memory location storing interleaved vectors with the same byte-widths as stored in widths
     * @param size size of the interleaved vectors stored at data
     */
    void set_data(char* data, uint64_t size) {
        data_vectors.clear();
        data_vectors.shrink_to_fit();
        
        size_vectors = size;
        capacity_vectors = size;
        
        set_bases(data);
    }

    /**
     * @brief serializes the interleaved vectors to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const {
        out.write((char*)&size_vectors,sizeof(uint64_t));
        out.write((char*)&width_entry,sizeof(uint64_t));

        if (num_vectors > 0) {
            out.write((char*)&widths[0],num_vectors*sizeof(uint64_t));
            out.write((char*)&masks[0],num_vectors*sizeof(uint_t));
        }

        if (size_vectors > 0) {
            write_to_file(out,(char*)&data_vectors[0],size_vectors*width_entry);
        }
    }

    /**
     * @brief loads the interleaved vectors from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        uint64_t old_size;

        in.read((char*)&old_size,sizeof(uint64_t));
        in.read((char*)&width_entry,sizeof(uint64_t));

        if (num_vectors > 0) {
            in.read((char*)&widths[0],num_vectors*sizeof(uint64_t));
            in.read((char*)&masks[0],num_vectors*sizeof(uint_t));
        }

        if (old_size > 0) {
            resize_no_init(old_size);
            read_from_file(in,(char*)&data_vectors[0],size_vectors*width_entry);
        }
    }

    std::ostream& operator>>(std::ostream& os) const {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is) {
        load(is);
        return is;
    }
};