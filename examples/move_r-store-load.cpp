#include <move_r/move_r.hpp>

int main()
{
    // build an index
    move_r<> index("This is a test string");

    // store an index in a file
    std::ofstream index_ofile("test_idx.move-r");
    index >> index_ofile;
    index_ofile.close();

    // load the same index into another move_r-object
    std::ifstream index_ifile("test_idx.move-r");
    move_r<> reloaded_index;
    reloaded_index << index_ifile;
    index_ifile.close();
}