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

template <int N, int K>
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
    var1 = double(var1) + double(var2);
    return var1;
}


template <typename T1, typename T2>
T1& operator -= (T1& var1, const T2& var2) {
    var1 = var1 - var2;
    return var1;
}

template <int w, int h>
class Size {
public:
    static constexpr int width = w;
    static constexpr int height = h;
    static std::string getName() {
        return "S(" + std::to_string(width) + "," + to_string(height) +")";
    }

};


template <typename PType, typename VType, typename VFlowType, int N, int M>
class Fluid{
public:
    static constexpr size_t T = 1'000'000;
    static constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    int rtN{};
    int rtM{};
    std::array<PType, 256> rho{};
    Fluid(const std::vector<std::string> &srcField) {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M + 1; j++) {
                field[i][j] = srcField[i][j];
            }
        }
    };
    Fluid(int rtN, int rtM): rtN(rtN), rtM(rtM), field(rtN, rtM + 1),
                             p(rtN, rtM), old_p(rtN, rtM), last_use(rtN, rtM), dirs(rtN, rtM), velocity(rtN, rtM), velocity_flow(rtN, rtM)  {};

    template <typename T, int planeN = N, int planeM = M>
    struct Plane {
        Plane() = default;
        Plane(int rtN, int rtM) : p(rtN, std::vector<T> (rtM)) {};
        using ptype = typename std::conditional_t<
                (N != 0), std::array<std::array<T, planeM>, planeN>,
        std::vector<std::vector<T>>
        >;

        ptype p;
    };

    Plane<char, N, M + 1> field;
    Plane<PType> p, old_p;
    Plane<int> last_use;
    Plane<int> dirs;
    template <typename T>
    struct VectorField {
        VectorField() = default;
        VectorField(int rtN, int rtM) : v(rtN, rtM){};
        Plane<array<T, deltas.size()>> v;

        T &add(int x, int y, int dx, int dy, T dv) {
            return get(x, y, dx, dy) += dv;
        }

        T &get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v.p[x][y][i];
        }
    };
    VectorField<VType> velocity{};
    VectorField<VFlowType>velocity_flow{};


    int UT = 0;


    mt19937 rnd{1337};

    tuple<VType, bool, pair<int, int>> propagate_flow(int x, int y, VType lim) {
        last_use.p[x][y] = UT - 1;
        VType ret = 0;
        for (auto [dx, dy]: deltas) {
            int nx = x + dx, ny = y + dy;
            if (field.p[nx][ny] != '#' && last_use.p[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = min(double(lim), double(cap - flow));
                if (last_use.p[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use.p[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, VFlowType(double(t)));
                    last_use.p[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                    return {t, prop && end != pair(x, y), end};
                }
            }
        }
        last_use.p[x][y] = UT;
        return {ret, 0, {0, 0}};
    }

    double random01() {
        std::uniform_real_distribution<double> dis(0.0, 1.0);
        return dis(rnd);
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy]: deltas) {
                int nx = x + dx, ny = y + dy;
                if (field.p[nx][ny] != '#' && last_use.p[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use.p[x][y] = UT;
        for (auto [dx, dy]: deltas) {
            int nx = x + dx, ny = y + dy;
            if (field.p[nx][ny] == '#' || last_use.p[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    VType move_prob(int x, int y) {
        VType sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field.p[nx][ny] == '#' || last_use.p[nx][ny] == UT) {
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
        PType cur_p;
        array<VType, deltas.size()> v;

        void swap_with(int x, int y, decltype(field) field, decltype(velocity) velocity, decltype(p) p) {
            swap(field.p[x][y], type);
            swap(p.p[x][y], cur_p);
            swap(velocity.v.p[x][y], v);
        }
    };

    bool propagate_move(int x, int y, bool is_first) {
        last_use.p[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<VType, deltas.size()> tres;
            VType sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field.p[nx][ny] == '#' || last_use.p[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0) {
                break;
            }

            VType p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0 && field.p[nx][ny] != '#' && last_use.p[nx][ny] < UT);

            ret = (last_use.p[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use.p[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field.p[nx][ny] != '#' && last_use.p[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(x, y, field, velocity, p);
                pp.swap_with(nx, ny, field, velocity, p);
                pp.swap_with(x, y, field, velocity, p);
            }
        }
        return ret;
    }

    int main() {
        rho[' '] = 0.01;
        rho['.'] = 1000;
        PType g = 0.1;

        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field.p[x][y] == '#')
                    continue;
                for (auto [dx, dy]: deltas) {
                    dirs.p[x][y] += (field.p[x + dx][y + dy] != '#');
                }
            }
        }

        for (size_t i = 0; i < T; ++i) {

            PType total_delta_p = 0;
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field.p[x][y] == '#')
                        continue;
                    if (field.p[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, VType(g));
                }
            }

            // Apply forces from p
            p = old_p;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field.p[x][y] == '#')
                        continue;
                    for (auto [dx, dy]: deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field.p[nx][ny] != '#' && old_p.p[nx][ny] < old_p.p[x][y]) {
                            auto delta_p = old_p.p[x][y] - old_p.p[nx][ny];
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field.p[nx][ny]] >= force) {
                                contr -= force / rho[(int) field.p[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field.p[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, force / rho[(int) field.p[x][y]]);
                            p.p[x][y] -= force / dirs.p[x][y];
                            total_delta_p -= force / dirs.p[x][y];
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
                        if (field.p[x][y] != '#' && last_use.p[x][y] != UT) {
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
                    if (field.p[x][y] == '#')
                        continue;
                    for (auto [dx, dy]: deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = double(new_v);
                            auto force = (old_v - new_v) * rho[(int) field.p[x][y]];
                            if (field.p[x][y] == '.')
                                force *= 0.8;
                            if (field.p[x + dx][y + dy] == '#') {
                                p.p[x][y] += force / dirs.p[x][y];
                                total_delta_p += force / dirs.p[x][y];
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
                    if (field.p[x][y] != '#' && last_use.p[x][y] != UT) {
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
                for (size_t y = 0; y < N; ++y) {
                    for(size_t x = 0; x < M + 1; ++x) {
                        cout << field.p[y][x];
                    }
                    cout << "\n";
                }
            }
        }
        return 0;
    }
};

#define SIZES S(36, 84),S(40,40)
#define S(w, h) Size<w, h>
#define TYPES FIXED(32,6),FAST_FIXED(30,7),DOUBLE
#define FIXED(a, b) Fixed<a,b>
#define FAST_FIXED(a, b) FastFixed<a,b>
#define DOUBLE Double
#define FLOAT Float
template <typename PType, typename VType, typename VFlowType, typename CurSize, typename... Sizes>
void getFluidSize(int width, int height, const std::vector<std::string>& field  ) {
    if(width == CurSize::width && height == CurSize::height) {
        auto fluid = Fluid<PType, VType, VFlowType, CurSize::width, CurSize::height>();
        fluid.main();
    }
    if constexpr (sizeof...(Sizes) > 0) {
        getFluidSize<PType, VType, VFlowType, Sizes...>(width, height, field);
    }
    ;

}
template <typename PType, typename VType, typename CurVFlowType, typename... Types>
auto getFluidVFlowType(const std::string& vFlowTypeName, int width, int height, const std::vector<std::string>& field) {
    if(vFlowTypeName == CurVFlowType::getName()) {
        getFluidSize<PType, VType, CurVFlowType, SIZES>(width, height, field);
    }
    if constexpr (sizeof...(Types) > 0) {
        getFluidVFlowType<PType, VType, Types...>(vFlowTypeName, width, height, field);
    }
    ;

}

template <typename PType, typename CurVType, typename... Types>
auto getFluidVType(const std::string& vTypeName, const std::string& vFlowTypeName, int width, int height, const std::vector<std::string>& field) {
    if(vTypeName == CurVType::getName()) {
        getFluidVFlowType<PType, CurVType, TYPES>(vFlowTypeName, width, height, field);
    }
    if constexpr (sizeof...(Types) > 0) {
        getFluidVType<PType, Types...>(vTypeName, vFlowTypeName, width, height, field);
    }
    ;

}

template <typename CurPType, typename... Types>
auto getFluid(const std::string& pTypeName, const std::string& vTypeName, const std::string& vFlowTypeName, int width, int height, const std::vector<std::string>& field) {
    if(pTypeName == CurPType::getName()) {
        getFluidVType<CurPType, TYPES>(vTypeName, vFlowTypeName, width, height, field);
    }
    if constexpr (sizeof...(Types) > 0) {
        getFluid<Types...>(pTypeName, vTypeName, vFlowTypeName, width, height, field);
    }
    ;

}

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
    std::vector<std::string > field;
    std::string line;

    while (std::getline(fin, line)) {
        field.push_back(line);
    }
    getFluid<TYPES>("DOUBLE", "DOUBLE", "DOUBLE", field.at(0).size() - 1, field.size(), field);

    return 0;
}
