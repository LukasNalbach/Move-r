#pragma once

#include <array>
#include <vector>
#include <iostream>

/**
 * @brief variable-width interleaved vectors (can store at most 16 interleaved vectors, widths are fixed to whole bytes)
 * @tparam uint_t unsigned integer type
 */
template <typename uint_t = uint32_t>
class interleaved_vectors {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    protected:
    uint64_t size = 0; // size of each stored vector
    uint8_t vecs = 0; // number of stored vectors
    uint64_t width_entry = 0; // sum of the widths of all vectors

    /**
     * @brief [0..(size+1)*vecs-1] vector storing the interleaved vectors
     */
    std::vector<uint8_t> data;

    /**
     * @brief [0..vecs-1] widths of the stored vectors; widths[i] = width of vector i
     */
    std::array<uint8_t,16> widths = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    /**
     * @brief [0..vecs-1] pointers to the first entries of each vector; bases[vec] = base
     *        of the first entry of the vector with index vec
     */
    std::array<uint8_t*,16> bases = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    /**
     * @brief [0..vecs-1] masks that are used to mask off data of other vector entries when
     *        accessing a vector
     */
    std::array<uint_t,16> masks_get = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    /**
     * @brief [0..vecs-1] masks that are used to mask off data of other vector entries when
     *        setting an entry in a vector
     */
    std::array<uint_t,16> masks_set = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    void set_bases() {
        bases[0] = &data[0];

        for (uint8_t i=1; i<vecs; i++) {
            bases[i] = bases[i-1] + widths[i-1];
        }
    }

    /**
     * @brief copies another interleaved_vectors object into this object
     * @param other another interleaved_vectors object
     */
    void copy_from_other(const interleaved_vectors& other) {
        size = other.size;
        vecs = other.vecs;
        width_entry = other.width_entry;
        data = other.data;
        widths = other.widths;
        masks_get = other.masks_get;
        masks_set = other.masks_set;

        set_bases();
    }

    /**
     * @brief moves another interleaved_vectors object into this object
     * @param other another interleaved_vectors object
     */
    void move_from_other(interleaved_vectors&& other) {
        size = std::move(other.size);
        vecs = std::move(other.vecs);
        width_entry = std::move(other.width_entry);
        data = std::move(other.data);
        widths = std::move(other.widths);
        bases = std::move(other.bases);
        masks_get = std::move(other.masks_get);
        masks_set = std::move(other.masks_set);
        set_bases();

        other.size = 0;
        other.vecs = 0;
        other.width_entry = 0;

        for (uint8_t i=0; i<16; i++) {
            other.bases[i] = NULL;
            other.widths[i] = 0;
            other.masks_get[i] = 0;
            other.masks_set[i] = 0;
        }
    }
    
    public:
    interleaved_vectors() = default;
    interleaved_vectors(interleaved_vectors&& other) {move_from_other(std::move(other));}
    interleaved_vectors(const interleaved_vectors& other) {copy_from_other(other);}
    interleaved_vectors& operator=(interleaved_vectors&& other) {move_from_other(std::move(other));return *this;}
    interleaved_vectors& operator=(const interleaved_vectors& other) {copy_from_other(other);return *this;}

    ~interleaved_vectors() {
        for (uint8_t i=0; i<vecs; i++) {
            bases[i] = NULL;
        }
    }

    /**
     * @brief Construct a new interleaved_vectors 
     * @param widths list containing the byte-widths of the interleaved arrays
     * @param size the number of entries to preallocate
     * @param initialize_memory controls whether the vectors should be initialized to 0
     */
    interleaved_vectors(std::vector<uint8_t> widths, uint64_t size = 0, bool initialize_memory = true) {
        this->vecs = widths.size();

        for (uint8_t i=0; i<vecs; i++) {
            this->widths[i] = widths[i];
            width_entry += widths[i];
            masks_get[i] = std::numeric_limits<uint_t>::max()>>(8*(sizeof(uint_t)-widths[i]));
            masks_set[i] = ~masks_get[i];
        }

        resize(size,initialize_memory);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() {
        return
            2*sizeof(uint64_t)+1+ // variables
            vecs+ // widths
            2*vecs*sizeof(uint_t)+ // masks
            (size+1)*width_entry; // data
    }

    /**
     * @brief returns number of interleaved vectors
     * @return the number of interleaved vectors
     */
    uint8_t num_vectors() {
        return vecs;
    }

    /**
     * @brief returns total width (number of bytes) per entry, that is the (sum of all widths)
     * @return number of bytes per entry
     */
    uint8_t bytes_per_entry() {
        return this->width_entry;
    }

    /**
     * @brief returns the width in bytes of the vector with index vec
     * @param vec vector index
     * @return its witdth in bytes
     */
    uint8_t width(uint8_t vec) {
        return widths[vec];
    }

    /**
     * @brief resizes all stored vectors to size
     * @param size size
     * @param initialize_memory controls whether the vectors should be initialized to 0
     */
    void resize(uint64_t size, bool initialize_memory = true) {
        this->size = size;

        if (initialize_memory) {
            data.resize((size+1)*width_entry);
        } else {
            (*reinterpret_cast<std::vector<no_init<char>>*>(&data)).resize((size+1)*width_entry);
        }
        
        set_bases();
    }

    /**
     * @brief (possibly re-)allocates the memory needed to resize all vectors to size
     * @param size size
     */
    void reserve(uint64_t size) {
        data.reserve((size+1)*width_entry);
        set_bases();
    }

    /**
     * @brief shrinks all vectors to size (possibly reallocates memory)
     */
    void shrink_to_fit() {
        data.shrink_to_fit();
        set_bases();
    }

    /**
     * @brief sets the i-th entry in the vector with index vec to v
     * @tparam vec vector index (0 <= vec < vecs)
     * @param i entry index (0 <= i < size)
     * @param value to store
     */
    template <uint8_t vec>
    inline void set(uint_t i, uint_t v) {
        static_assert(vec < 16);
        *reinterpret_cast<uint_t*>(bases[vec] + i * width_entry) =
        (*reinterpret_cast<uint_t*>(bases[vec] + i * width_entry) & masks_set[vec]) | v;
    }

    /**
     * @brief plainly stores the value v (of type T) at the memory location of the i-th entry in the vector
     * with index vec (if sizeof(T) > widths[vec], this is unsafe, since subsequnet entries are overwritten)
     * @tparam vec vector index (0 <= vec < vecs)
     * @tparam T type to treat the i-th entry in the vector with index vec as
     * @param i entry index (0 <= i < size)
     * @param value to store
     */
    template <uint8_t vec, typename T>
    inline void set_unsafe(uint_t i, T v) {
        static_assert(vec < 16);
        *reinterpret_cast<T*>(bases[vec] + i * width_entry) = v;
    }

    /**
     * @brief returns the i-th entry in the vector with index vec
     * @tparam vec vector index (0 <= vec < vecs)
     * @param i entry index (0 <= i < size)
     * @return value
     */
    template <uint8_t vec>
    inline uint_t get(uint_t i) {
        static_assert(vec < 16);
        return *reinterpret_cast<uint_t*>(bases[vec] + i * width_entry) & masks_get[vec];
    }

    /**
     * @brief interprets the memory location of the i-th entry in the vector with index vec as an object of type
     * T and returns it (if sizeof(T) > widths[vec], the returned value contains data from subsequent entries)
     * @tparam vec vector index (0 <= vec < vecs)
     * @tparam uint_t type to treat the i-th entry in the vector with index vec as
     * @param i entry index (0 <= i < size)
     * @return value
     */
    template <uint8_t vec, typename T>
    inline uint_t get_unsafe(uint_t i) {
        static_assert(vec < 16);
        return *reinterpret_cast<T*>(bases[vec] + i * width_entry);
    }

    /**
     * @brief serializes the interleaved vectors to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) {
        out.write((char*)&size,sizeof(uint64_t));
        out.write((char*)&vecs,1);
        out.write((char*)&width_entry,sizeof(uint64_t));

        if (vecs > 0) {
            out.write((char*)&widths[0],vecs);
            out.write((char*)&masks_get[0],vecs*sizeof(uint_t));
            out.write((char*)&masks_set[0],vecs*sizeof(uint_t));
        }

        if (size > 0) {
            out.write((char*)&data[0],size*width_entry);
        }
    }

    /**
     * @brief loads the interleaved vectors from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&vecs,1);
        in.read((char*)&width_entry,sizeof(uint64_t));

        if (vecs > 0) {
            in.read((char*)&widths[0],vecs);
            in.read((char*)&masks_get[0],vecs*sizeof(uint_t));
            in.read((char*)&masks_set[0],vecs*sizeof(uint_t));
        }

        if (size > 0) {
            resize(size);
            in.read((char*)&data[0],size*width_entry);
        }
    }
};