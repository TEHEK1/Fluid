#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <any>
#include <cassert>
#include <random>
#include <string>
#include <fstream>
#include <memory>

using namespace std;

constexpr size_t N = 36, M = 84;
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1'000'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

class Type {
public:
};

class Float : public Type{
public:
    Float(float value) : m_value(value) {}
    Float() : Float(0) {};
    operator double () const{
        return m_value;
    }
    static std::string getName() {
        return "FLOAT";
    }
private:
    float m_value;
};

class Double : public Type{
public:
    Double(double value) : m_value(value) {}
    Double() : Double(0) {};
    operator double () const {
        return m_value;
    }
    static std::string getName() {
        return "DOUBLE";
    }
private:
    double m_value;
};

template <int N = 32, int K = 16>
class Fixed : public Type {
public:
    using ValueType = typename std::conditional_t<
            (N == 8), int8_t,
            typename std::conditional_t<
                    (N == 16), int16_t,
                    typename std::conditional_t<
                            (N == 32), int32_t,
                            typename std::conditional_t<
                                    (N == 64), int64_t,
                                    nullptr_t
                            >
                    >
            >
    >;
    static_assert(!std::is_same<ValueType, nullptr_t>::value);
    Fixed(double value) : m_value(value * (1 << K)) {
    }
    Fixed() : Fixed(0) {};
    operator double () const {
        return static_cast<double>(m_value) / (1 << K);
    }
    static std::string getName() {
        return "FIXED(" + std::to_string(N) + "," + to_string(K) +")";
    }
    ValueType m_value;
};

template <int N, int K>
class FastFixed : public Type {
public:
    static_assert(N <= 64, "FastFixed requires N <= 64");
    using ValueType = typename std::conditional_t<
            (N <= 8), int_fast8_t,
            typename std::conditional_t<
                    (N <= 16), int_fast16_t,
                    typename std::conditional_t<
                            (N <= 32), int_fast32_t,
                            typename std::conditional_t<
                                    (N <= 64), int_fast64_t,
                                    nullptr_t
                            >
                    >
            >
    >;
    static_assert(!std::is_same<ValueType, nullptr_t>::value);
    FastFixed(double value) : m_value(value * (1 << K)) {
    }
    FastFixed() : FastFixed(0) {};
    operator double () const {
        return static_cast<double>(m_value) / (1 << K);
    }
    static std::string getName() {
        return "FAST_FIXED(" + std::to_string(N) + "," + to_string(K) +")";
    }
    ValueType m_value;
};

template <typename T1, typename T2>
T1& operator += (T1& var1, const T2& var2) {
    var1 = var1 + var2;
    return var1;
}


template <typename T1, typename T2>
T1& operator -= (T1& var1, const T2& var2) {
    var1 = var1 - var2;
    return var1;
}

class Fluid{
public:
    Fluid(const std::vector<std::string> &srcField) {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                field.p.at(i).at(j) = srcField.at(i).at(j);
            }
            field.p.at(i).at(M) = '\0';
        }
    };

// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };
    template <typename T, int planeN = N, int planeM = M>
    struct Plane {
        Plane() = default;
        Plane(int rtN, int rtM) : p(rtN, std::vector<T> (rtM)) {};
        Plane(std::array<std::array<char, planeM>, planeN> srcfield) : p(srcfield) {
        }
        using ptype = typename std::conditional_t<
                (N != 0), std::array<std::array<T, planeM>, planeN>,
                std::vector<std::vector<T>>
        >;
        ptype p;
    };

    Plane<char, N, M+1> field{};




    std::array<Fixed<>, 256> rho;

    Plane<Fixed<>> p, old_p;


    struct VectorField {
        Plane<array<Fixed<>, deltas.size()>> v;
        Fixed<> &add(int x, int y, int dx, int dy, Fixed<> dv) {
            return get(x, y, dx, dy) += dv;
        }

        Fixed<> &get(int x, int y, int dx, int dy) {
            if(dx == -1 && dy == 0) { // Changed ranges to if
                return v.p.at(x).at(y).at(0);
            }
            else if(dx == 1 && dy == 0) {
                return v.p.at(x).at(y).at(1);
            }
            else if(dx == 0 && dy == -1) {
                return v.p.at(x).at(y).at(2);
            }
            else if(dx == 0 && dy == -1) {
                return v.p.at(x).at(y).at(2);
            }
            else {
                assert(false);
            }
        }
    };

    VectorField velocity{}, velocity_flow{};
    Plane<int> last_use;
    int UT = 0;


    mt19937 rnd{1337};

    tuple<Fixed<>, bool, pair<int, int>> propagate_flow(int x, int y, Fixed<> lim) {
        last_use.p.at(x).at(y) = UT - 1;
        Fixed<> ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field.p.at(nx).at(ny) != '#' && last_use.p.at(nx).at(ny)< UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = min(lim, Fixed<>(cap - flow));
                if (last_use.p.at(nx).at(ny) == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use.p.at(x).at(y) = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use.p.at(x).at(y) = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                    return {t, prop && end != pair(x, y), end};
                }
            }
        }
        last_use.p.at(x).at(y) = UT;
        return {ret, 0, {0, 0}};
    }

    double random01() {
        std::uniform_real_distribution<double> dis(0.0, 1.0);
        return dis(rnd);
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field.p.at(nx).at(ny)!= '#' && last_use.p.at(nx).at(ny) < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use.p.at(x).at(y) = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field.p.at(nx).at(ny) == '#' || last_use.p.at(nx).at(ny) == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    Fixed<> move_prob(int x, int y) {
        Fixed<> sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field.p.at(nx).at(ny) == '#' || last_use.p.at(nx).at(ny) == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    struct ParticleParams {
        char type;
        Fixed<> cur_p;
        array<Fixed<>, deltas.size()> v;

        void swap_with(int x, int y, Plane<char, N, M + 1>& field, Plane<Fixed<>>& p, VectorField& velocity) {
            swap(field.p.at(x).at(y), type);
            swap(p.p.at(x).at(y), cur_p);
            swap(velocity.v.p.at(x).at(y), v);
        }
    };

    bool propagate_move(int x, int y, bool is_first) {
        last_use.p.at(x).at(y) = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<Fixed<>, deltas.size()> tres;
            Fixed<> sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field.p.at(nx).at(ny) == '#' || last_use.p.at(nx).at(ny) == UT) {
                    tres.at(i) = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres.at(i) = sum;
                    continue;
                }
                sum += v;
                tres.at(i) = sum;
            }

            if (sum == 0) {
                break;
            }

            Fixed<> p = random01() * sum;
            size_t d = tres.at(0) > p ? 0 : tres.at(1) > p ? 1 : tres.at(2) > p ? 2 :tres.at(3) > p ? 3 : 4; // Changed ranges to if

            auto [dx, dy] = deltas.at(d);
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0 && field.p.at(nx).at(ny) != '#' && last_use.p.at(nx).at(ny) < UT);

            ret = (last_use.p.at(nx).at(ny) == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use.p.at(nx).at(ny) = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas.at(i);
            int nx = x + dx, ny = y + dy;
            if (field.p.at(nx).at(ny) != '#' && last_use.p.at(nx).at(ny) < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(x, y, field, p, velocity);
                pp.swap_with(nx, ny, field, p, velocity);
                pp.swap_with(x, y, field, p, velocity);
            }
        }
        return ret;
    }

    Plane<int> dirs;

    int main() {
        rho[' '] = 0.01;
        rho['.'] = 1000;
        Fixed<> g = 0.1;

        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field.p.at(x).at(y) == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs.p.at(x).at(y) += (field.p.at(x + dx).at(y + dy) != '#');
                }
            }
        }

        for (size_t i = 0; i < T; ++i) {

            Fixed<> total_delta_p = 0;
            std::array<size_t, N> xRange;
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field.p.at(x).at(y) == '#')
                        continue;
                    if (field.p.at(x + 1).at(y) != '#')
                        velocity.add(x, y, 1, 0, g);
                }
            }

            // Apply forces from p
            old_p = p;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field.p.at(x).at(y) == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field.p.at(nx).at(ny) != '#' && old_p.p.at(nx).at(ny) < old_p.p.at(x).at(y)) {
                            auto delta_p = old_p.p.at(x).at(y) - old_p.p[nx][ny];
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field.p[nx][ny]] >= force) {
                                contr -= force / rho[(int) field.p[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field.p[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, force / rho[(int) field.p[x][y]]);
                            p.p.at(x).at(y) -= force / dirs.p.at(x).at(y);
                            total_delta_p -= force / dirs.p.at(x).at(y);
                        }
                    }
                }
            }

            // Make flow from velocities
            velocity_flow = {};
            bool prop = false;
            do {
                UT += 2;
                prop = 0;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field.p.at(x).at(y) != '#' && last_use.p.at(x).at(y) != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = 1;
                            }
                        }
                    }
                }
            } while (prop);

            // Recalculate p with kinetic energy
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field.p.at(x).at(y) == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = (old_v - new_v) * rho[(int) field.p[x][y]];
                            if (field.p.at(x).at(y) == '.')
                                force *= 0.8;
                            if (field.p[x + dx][y + dy] == '#') {
                                p.p.at(x).at(y) += force / dirs.p.at(x).at(y);
                                total_delta_p += force / dirs.p.at(x).at(y);
                            } else {
                                p.p[x + dx][y + dy] += force / dirs.p[x + dx][y + dy];
                                total_delta_p += force / dirs.p[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field.p.at(x).at(y) != '#' && last_use.p.at(x).at(y) != UT) {
                        if (random01() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < N; ++x) {
                    for(size_t y = 0; y < field.p[x].size(); y++) {
                        if(!x) {
                            break;
                        }
                        cout << field.p.at(x).at(y);
                    }
                    cout <<  "\n";
                }
            }
        }
    }
};

int main(int argc, char* argv[]) {
    std::string pTypeName;
    std::string vTypeName;
    std::string vFlowTypeName;
    std::string fileName;
    std::vector< std::vector<char> > startFame;
    for(int i = 1; i < argc; i++) {
        std::string cur_arg(argv[i]);
        if(cur_arg.starts_with("--p-type=")) {
            pTypeName = cur_arg.substr(cur_arg.find('='));
        }
        else if(cur_arg.starts_with("--v-type=")) {
            vTypeName = cur_arg.substr(cur_arg.find('='));
        }
        else if(cur_arg.starts_with("--v-flow-type=")) {
            vFlowTypeName = cur_arg.substr(cur_arg.find('='));
        }
        else {
            fileName = cur_arg;
        }
    }
    assert(!fileName.empty());

    std::ifstream fin{fileName};
    if(!fin) {
        std::cout << "Failed to open " + fileName;
        exit(EXIT_FAILURE);
    }
    std::vector<std::string > field;
    std::string line;
    while (std::getline(fin, line)) {
        field.push_back(line);
    }
    assert(field.size());

    Fluid fluid(field);
    fluid.main();
}