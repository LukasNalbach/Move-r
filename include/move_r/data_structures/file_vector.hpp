#pragma once

/**
 * @brief interface, to read from and write to a file like a vector
 * @tparam val_t value type
 * @tparam pos_t position type
 */
template <typename val_t, typename pos_t>
class file_vector {
    protected:
    std::ifstream* ifile = NULL;
    std::ofstream* ofile = NULL;

    uint64_t width = 0;
    uint64_t pos = 0;

    public:
    file_vector() = default;

    file_vector(std::ifstream& ifile, uint8_t width) {
        this->ifile = &ifile;
        this->width = width;
    }
    
    file_vector(std::ofstream& ofile, uint8_t width) {
        this->ofile = &ofile;
        this->width = width;
    }

    inline pos_t operator[](pos_t i) {
        ifile->seekg((uint64_t{i}*width)-pos,std::ios_base::cur);
        pos_t val = 0;
        ifile->read((char*)&val,width);
        pos = (uint64_t{i}+1)*width;
        return val;
    }

    inline void write(pos_t i, val_t val) {
        ofile->seekp((uint64_t{i}*width)-pos,std::ios_base::cur);
        ofile->write((char*)&val,width);
        pos = (uint64_t{i}+1)*width;
    }
};