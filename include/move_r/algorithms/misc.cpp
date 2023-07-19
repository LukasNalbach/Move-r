template <typename uint_t>
template <typename output_t, bool output_reversed>
void move_r<uint_t>::retrieve_range(
    void(move_r<uint_t>::*retrieve_method)(const std::function<void(uint_t,output_t)>&,uint_t,uint_t,uint16_t),
    std::ofstream& out, uint_t l, uint_t r, uint16_t num_threads, int64_t max_bytes_alloc
) {
    int64_t max_writebuffer_size = max_bytes_alloc != -1 ? max_bytes_alloc/(num_threads*sizeof(output_t)) : std::max((int64_t)1,(int64_t)(r-l+1)/(num_threads*500));

    if (num_threads == 1 && !output_reversed) {
        std::vector<output_t> write_buffer;
        (*reinterpret_cast<std::vector<no_init<output_t>>*>(&write_buffer)).resize(max_writebuffer_size);
        int64_t pos_cur;

        (this->*retrieve_method)([&out,&pos_cur,&write_buffer,&max_writebuffer_size](uint_t, output_t v){
            write_buffer[pos_cur++] = v;
            
            if (pos_cur == max_writebuffer_size) {
                out.write((char*)&write_buffer[0],max_writebuffer_size*sizeof(output_t));
                pos_cur = 0;
            }
        },l,r,num_threads);
    } else {
        std::vector<std::ofstream> files;
        std::string file_prefix = "move-r_" + random_alphanumeric_string(10) + "_";
        std::vector<std::vector<output_t>> write_buffer(num_threads);
        
        for (uint16_t i=0; i<num_threads; i++) {
            files.emplace_back(std::ofstream(file_prefix + std::to_string(i)));
        }

        #pragma omp parallel for num_threads(num_threads)
        for (uint_t i=0; i<num_threads; i++) {
            (*reinterpret_cast<std::vector<no_init<output_t>>*>(&write_buffer[i])).resize(max_writebuffer_size);
        }
        
        std::vector<int64_t> pos_cur(num_threads,0);
        
        (this->*retrieve_method)([&files,&pos_cur,&write_buffer,&max_writebuffer_size](uint_t, output_t v){
            write_buffer[omp_get_thread_num()][pos_cur[omp_get_thread_num()]++] = v;
            
            if (pos_cur[omp_get_thread_num()] == max_writebuffer_size) {
                files[omp_get_thread_num()].write((char*)&write_buffer[omp_get_thread_num()][0],max_writebuffer_size*sizeof(output_t));
                pos_cur[omp_get_thread_num()] = 0;
            }
        },l,r,num_threads);

        #pragma omp parallel for num_threads(num_threads)
        for (uint16_t i=0; i<num_threads; i++) {
            if (pos_cur[i] > 0) {
                files[i].write((char*)&write_buffer[i][0],pos_cur[i]*sizeof(output_t));
            }
        }

        write_buffer.clear();
        write_buffer.shrink_to_fit();

        pos_cur.clear();
        pos_cur.shrink_to_fit();

        for (uint16_t i=0; i<num_threads; i++) {
            files[i].close();
        }

        files.clear();

        std::vector<output_t> temp_buffer;
        int64_t max_tempbuffer_size = num_threads*max_writebuffer_size;
        (*reinterpret_cast<std::vector<no_init<output_t>>*>(&temp_buffer)).resize(max_tempbuffer_size);
        int64_t elems_left;
        int64_t elems_to_read;

        for (uint16_t i=0; i<num_threads; i++) {
            std::ifstream file(file_prefix + std::to_string(i));
            file.seekg(0,std::ios::end);
            elems_left = file.tellg()/sizeof(output_t);
            if constexpr (!output_reversed) file.seekg(0,std::ios::beg);

            while (elems_left > 0) {
                elems_to_read = std::min(elems_left,max_tempbuffer_size);
                if constexpr (output_reversed) file.seekg(file.tellg()-(int64_t)(elems_to_read*sizeof(output_t)));
                file.read((char*)&temp_buffer[0],elems_to_read*sizeof(output_t));
                if constexpr (output_reversed) std::reverse(&temp_buffer[0],&temp_buffer[elems_to_read]);
                out.write((char*)&temp_buffer[0],elems_to_read*sizeof(output_t));
                if constexpr (output_reversed) file.seekg(file.tellg()-(int64_t)(elems_to_read*sizeof(output_t)));
                elems_left -= elems_to_read;
            }

            file.close();
            std::filesystem::remove(file_prefix + std::to_string(i));
        }
    }
}