#ifndef UTIL_H
#define UTIL_H


template <typename T>
void save_basic_type(std::ostream& os, const T& val, int version = 1)
{
    if (version != 1)
    {
        std::cerr << "version not implemented!\n";
        exit(EXIT_FAILURE);
    }
    uint64_t tmp = uint64_t(val);
    os.write((char*)&tmp, 8);
}

template <typename T>
void load_basic_type(std::istream& is, T& val)
{
    uint64_t tmp;
    is.read((char*)&tmp, 8);
    if (not is.good())
    {
        std::cerr << "error reading basic type\n";
        exit(EXIT_FAILURE);
    }
    val = T(tmp);
}

// version 1:
//   version, size(), raw array
template <typename T>
void save_shallow_vector(std::ostream& os, const std::vector<T>& v, int version = 1)
{
    if (version != 1)
    {
        std::cerr << "version not implemented!\n";
        exit(EXIT_FAILURE);
    }

    save_basic_type<uint64_t>(os, version, 1);
    save_basic_type<uint64_t>(os, v.size(), 1);
    os.write((char*)&v[0], v.size() * sizeof(T));
}

template <typename T>
void load_shallow_vector(std::istream& is, std::vector<T>& v)
{
    int version;
    load_basic_type(is, version);
    if (version != 1)
    {
        std::cerr << "version not implemented!\n";
        exit(EXIT_FAILURE);
    }
    size_t n;
    load_basic_type(is, n);
    v.resize(n);
    is.read((char*)&v[0], v.size() * sizeof(T));
    if (not is.good())
    {
        std::cerr << "error reading shallow vector\n";
        exit(EXIT_FAILURE);
    }
}


#endif
