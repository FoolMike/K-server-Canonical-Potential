#include <cstdio>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <queue>
#include <cmath>
#include <complex>
#include <bitset>
using namespace std;

const double inf = 1e9, eps = 1e-9;
const int LENGTH = 6;

bool flag;
double range;

template<class T>
vector<T> Config(T x, T y, T z)
{
    vector<T> config;
    config.push_back(x);
    config.push_back(y);
    config.push_back(z);
    return config;
}
template<class T>
pair<vector<T>, T> SupportPoint(T x, T y, T z, double w)
{
    vector<T> config;
    config.push_back(x);
    config.push_back(y);
    config.push_back(z);
    return make_pair(config, w);
}

template<int K, class Point>
class WorkFunctionSimulator{
public:
    WorkFunctionSimulator(double (*distance)(Point, Point)): distance(distance), ratio1(0), ratio2(-1), ratio3(-0.5), ratio4(0), ratio5(0) {}

    typedef vector<Point> Configuration;
    typedef vector<pair<Configuration, double> > Support;

    double (*distance)(Point, Point);
    Point (*antipode)(Point);
    void (*print)(Configuration);
    // int vis[K], match[K], stamp;
    // double lx[K], ly[K], delta;
    vector<Point> enumerator;

    bool in(Configuration C, Point P)
    {
        for (Point p: C)
            if (distance(p, P) < eps)
                return true;
        return false;
    }

    // bool find(int x)
    // {
    //     vis[x] = stamp;
    //     for (int i = 0; i < K; ++i)
    //     {
    //         if (abs(lx[x] + ly[i] - f[x][i]) < eps)
    //         {
    //             if (match[i] == -1 || (vis[match[i]] != stamp && find(match[i])))
    //             {
    //                 match[i] = x;
    //                 return true;
    //             }
    //         }
    //         if (f[x][i] - lx[x] - ly[i] > eps && (~match[i] || vis[match[i]] != stamp))
    //         {
    //             delta = min(delta, f[x][i] - lx[x] - ly[i]);
    //         }
    //     }
    //     return false;
    // }

    double distanceC(Configuration A, Configuration B) const
    {
        if (K <= 6) return _distanceC(A, B);
        // for (int i = 0; i < K; ++i)
        // {
        //     lx[i] = inf;
        //     ly[i] = 0;
        //     for (int j = 0; j < K; ++j)
        //     {
        //         f[i][j] = distance(A[i], B[j]);
        //         lx[i] = min(lx[i], f[i][j]);
        //     }
        // }
        // memset(match, -1, sizeof(match));

        // for (int i = 0; i < K; ++i)
        // {
        //     while (1)
        //     {
        //         ++stamp;
        //         delta = inf;
        //         if (find(i)) break;
        //         for (int i = 0; i < K; ++i) 
        //         {
        //             if (vis[i] == stamp) lx[i] += delta;
        //             if (vis[match[i]] == stamp) ly[i] -= delta;
        //         }
        //     }
        // }

        // double cost = 0;
        // for (int i = 0; i < K; ++i)
        // {
        //     //assert(~match[i]);
        //     cost += f[match[i]][i];
        // }
        // return cost;
    }

    double _distanceC(Configuration A, Configuration B) const
    {
        static int perm[K];
        for (int i = 0; i < K; ++i) perm[i] = i;
        double ans = inf;
        static double f[K][K];
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j)
                f[i][j] = distance(A[i], B[j]);
        do{
            double cost = 0;
            for (int i = 0; i < K; ++i) cost += f[i][perm[i]];
            ans = min(ans, cost);
        }while (next_permutation(perm, perm + K));
        return ans;
    }

    Configuration C;
    vector<Point> requests;
    Support supp;
    double Cost, ExtendedCost;

    double workFunction(Configuration A) const
    {
        double cost = inf;
        for (pair<Configuration, double> C: supp)
            cost = min(cost, C.second + distanceC(A, C.first));
        return cost;
    }

    double workFunction(Configuration A, const Support& supp) const
    {
        double cost = inf;
        for (pair<Configuration, double> C: supp)
            cost = min(cost, C.second + distanceC(A, C.first));
        return cost;
    }

    Support append(Support supp, Point r)
    {
        Support ans, all;
        for (int i = 0; i < supp.size(); ++i)
        {
            Configuration A = supp[i].first;
            for (int j = 0; j < K; ++j)
            {
                swap(r, A[j]);
                all.push_back(make_pair(A, supp[i].second + distance(A[j], r)));
                swap(r, A[j]);
            }
        }
        for (int i = 0; i < all.size(); ++i)
        {
            bool ins = true;
            for (int j = 0; j < all.size(); ++j)
                if (i != j)
                {
                    double dis = distanceC(all[i].first, all[j].first);
                    if (all[i].second > all[j].second && all[i].second >= all[j].second + dis - eps)
                    {
                        ins = false;
                        break;
                    }
                }
            if (ins && workFunction(all[i].first, ans) > all[i].second) ans.push_back(all[i]);
        }
        return ans;
    }

    void init(Configuration A)
    {
        C = A;
        Cost = ExtendedCost = 0;
        requests.clear();
        supp.clear();
        supp.push_back(make_pair(A, 0));
    }

    void append(Point r)
    {
        Configuration Nxt, Ext;
        for (int i = 0; i < K; ++i) Ext.push_back(antipode(r));
        double potential = inf, wf = workFunction(Ext);
        // for (int i = 0; i < K; ++i)
        // {
        //     swap(C[i], r);
        //     double cost = workFunction(C, supp) + distance(C[i], r);
        //     if (cost < potential)
        //     {
        //         potential = cost;
        //         Nxt = C;
        //     }
        //     swap(C[i], r);
        // }
        // Cost += distanceC(Nxt, C);
        // C = Nxt;
        requests.push_back(r);
        supp = append(supp, r);
        double extra_ext = workFunction(Ext) - wf;
        ExtendedCost += extra_ext;
        /*printf("cost = %lf\n", distanceC(Nxt, C));
        printf("extended cost = %lf\n", extra_ext);
        printf("New Configuration is\n");
        print(C);
        puts("\n");*/
    }

    void PrintSupport() const
    {
        for (pair<Configuration, double> C: supp)
        {
            print(C.first);
            printf(" w = %lf\n", C.second);
        }
        puts("");
    }

    double optimal() const
    {
        double ans = inf;
        for (pair<Configuration, double> C: supp)
        {
            ans = min(ans, C.second);
        }
        return ans;
    }

    double clique(Configuration C) const
    {
        double ans = 0;
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < i; ++j)
                ans += distance(C[i], C[j]);
        return ans;
    }

    double ratio1, ratio2, ratio3, ratio4, ratio5;

    double Potential(Point u, Point a, Point y, Point p, Point b, Point x, Point q) const
    {
        return min(Potential_3Circle_Canonical(u, a, y, p, b, x, q), Potential_3Tree_Canonical(u, a, y, p, b, x, q));
    }

    double Potential_New(Point u, Point a, Point y, Point p, Point b, Point x, Point q, int Edge[LENGTH]) const
    {
        Point P[10] = {a, y, p, b, x, q, u};
        Point baru = antipode(u);
        Configuration C0 = Config(baru, baru, baru), C1 = Config(u, a, y), C2 = Config(u, p, b), C3 = Config(u, x, q);
        double value = workFunction(C0) + workFunction(C1) + workFunction(C2) + workFunction(C3);
        for (int i = 0; i < LENGTH; ++i)
            for (int j = 0; j < i; ++j)
                if ((Edge[i] >> j) & 1) value -= distance(P[i], P[j]);
        return value;
    }
    double Potential_New(Point u, int Edge[LENGTH]) const
    {
        double value = inf;
        for (Point a: enumerator)
            for (Point y: enumerator)
                for (Point p: enumerator)
                    for (Point b: enumerator)
                        for (Point x: enumerator)
                            for (Point q: enumerator)
                                value = min(value, Potential_New(u, a, y, p, b, x, q, Edge));
        return value;
    }
    double Potential_New(int Edge[LENGTH]) const
    {
        double value = inf;
        for (Point u: enumerator)
            value = min(value, Potential_New(u, Edge));
        return value;
    }

    double Potential_3Circle_Canonical(Point u, Point a, Point y, Point p, Point b, Point x, Point q) const
    {
        //assert(K == 3);
        Point baru = antipode(u);
        Configuration C0 = Config(baru, baru, baru), C1 = Config(u, a, y), C2 = Config(u, p, b), C3 = Config(u, x, q);
        double value = workFunction(C0) + workFunction(C1) + workFunction(C2) + workFunction(C3);
        value -= distance(a, b) + distance(p, q) + distance(x, y);
        value -= distance(a, p) + distance(p, x) + distance(x, a);
        //value -= distance(y, b) + distance(b, q) + distance(q, y);
        return value;
    }

    double Potential_3Tree_Canonical(Point u, Point a, Point y, Point p, Point b, Point x, Point q) const
    {
        //assert(K == 3);
        Point baru = antipode(u);
        Configuration C0 = Config(baru, baru, baru), C1 = Config(u, a, y), C2 = Config(u, p, b), C3 = Config(u, x, q);
        double value = workFunction(C0) + workFunction(C1) + workFunction(C2) + workFunction(C3);
        value -= distance(a, b);
        value -= distance(a, p) + distance(x, a);
        value -= distance(b, q) + distance(y, b);
        //value -= distance(a, p) + distance(p, x) + distance(x, a);
        //value -= distance(b, q) + distance(q, y) + distance(y, b);
        return value;
    }

    double Potential_3Tree_Alternative(Point u, Point x, Point y) const
    {
        //assert(K == 3);
        Point baru = antipode(u);
        Point barx = antipode(x);
        Point bary = antipode(y);
        Configuration C0 = Config(baru, baru, baru), C1 = Config(u, x, y), C2 = Config(u, x, bary), C3 = Config(u, barx, y);
        double value = workFunction(C0) + workFunction(C1) + workFunction(C2) + workFunction(C3) + distance(barx, y);
        return value;
    }

    double Potential_Another(Point u, Point a, Point y, Point p, Point b, Point x, Point q) const
    {
        //assert(K == 3);
        Point baru = antipode(u);
        Configuration C0 = Config(baru, baru, baru), C1 = Config(u, a, y), C2 = Config(u, p, b), C3 = Config(u, x, q);
        double value = workFunction(C0) + workFunction(C1) + workFunction(C2) + workFunction(C3);
        value -= distance(a, b) + distance(b, q) + distance(q, p) + distance(p, x) + distance(x, y);
        return value;
    }

    double Potential_3Circle(Point u, Point x, Point y, Point z) const
    {
        //assert(K == 3);
        Point barw = antipode(w);
        Point barx = antipode(x);
        Point bary = antipode(y);
        Point barz = antipode(z);
        double value = 0;
        value = workFunction(Config(baru, baru, baru))
                + workFunction(Config(u, x, bary))
                + workFunction(Config(u, y, barz))
                + workFunction(Config(u, z, barx));
        value -= distance(x, y);
        value -= distance(y, z);
        value -= distance(z, x);
        return value;
    }

    double Potential_Zhiyi(Point x[K + 1][K + 1]) const
    {
        Configuration C[K + 1];
        for (int i = 0; i <= K; ++i)
            for (int j = i + 1; j <= K; ++j)
                x[j][i] = antipode(x[i][j]);
        for (int i = 0; i <= K; ++i)
            for (int j = 0; j <= K; ++j)
                if (i != j) C[i].push_back(x[i][j]);
        double value = 0;
        for (int i = 0; i <= K; ++i) value += workFunction(C[i]) + clique(C[i]);
        return value;
    }

    double Potential_CK21(Point x, Point y, Point z) const
    {
        //assert(K == 3);
        Point barx = antipode(x);
        Point bary = antipode(y);
        Point barz = antipode(z);
        return workFunction(Config(x, y, z)) + workFunction(Config(barx, y, z)) + 
               workFunction(Config(bary, bary, z)) + workFunction(Config(barz, barz, barz));
    }

    Point Resolve(Configuration C) const
    {
        Point r = *requests.rbegin();
        print(C);
        printf("\n");
        double w = workFunction(C);
        for (int i = 0; i < K; ++i)
        {
            swap(C[i], r);
            if (workFunction(C) + distance(C[i], r) < w + eps)
            {
                cout << "Resolves from " << r << endl;
                //return r;
            }
            swap(C[i], r);
        }
        return r;
    }

    Configuration ResolveC(Configuration C) const
    {
        double w = workFunction(C);
        print(C);
        for (auto X: supp)
        {
            if (distanceC(C, X.first) + X.second == w)
            {
                cout << " Resolve from ";
                print(X.first);
                cout << endl;
                return X.first;
            }
        }
        return Configuration();
    }

    bool valid() const
    {
        Point r = *requests.rbegin();
        for (auto X: supp)
        {
            Point a = X.first[1], b = X.first[2];
            for (auto Y: supp)
            {
                Point c = Y.first[1], d = Y.first[2];
                double total = X.second + Y.second;
                if (workFunction(Config(r, a, c)) + workFunction(Config(r, b, d)) < total + eps) continue;
                if (workFunction(Config(r, a, d)) + workFunction(Config(r, b, c)) < total + eps) continue;
                return false;
            }
        }
        return true;
    }
};

double LineMetricDistance(double x, double y)
{
    return abs(x - y);
}
double CycleLength = 6;
double CycleMetricDistance(double x, double y)
{
    if (x < y) swap(x, y);
    return min(x - y, CycleLength - x + y);
}
double CycleMetricAntipode(double x)
{
    return x >= CycleLength / 2 ? x - CycleLength / 2 : x + CycleLength / 2;
}
template<class T>
void print(vector<T> C)
{
    cout << "(";
    for (unsigned int i = 0; i < C.size(); ++i)
    {
        cout << C[i];
        if (i + 1 != C.size()) cout << ' ';
    }
    cout << ")";
}

const int NumLeaves = 100;
double weight[NumLeaves];
double WeightedStarMetricDistance(int x, int y)
{
    return x == y ? 0 : weight[x] + weight[y];
}

double EuclideanMetricDistance(complex<double> x, complex<double> y)
{
    return norm(x - y);
}

double rand01()
{
    return 1.0 * rand() / (1ll << 31);
}

WorkFunctionSimulator<3, double> LineSimulator(LineMetricDistance);
WorkFunctionSimulator<3, double> CycleSimulator(CycleMetricDistance), _CycleSimulator(CycleMetricDistance), CycleSimulatorR(CycleMetricDistance);
WorkFunctionSimulator<3, int> WeightedStarSimulator(WeightedStarMetricDistance);
WorkFunctionSimulator<3, complex<double> > EuclideanPlaneSimulator(EuclideanMetricDistance);

const int K = 3;
vector<double> conf;
double r;
char command[10];

WorkFunctionSimulator<3, double> *Simulator = &CycleSimulator;

int n[2], m[2], type;
double dis[2][1010][1010];

double GraphMetricDistance(int x, int y)
{
    return dis[type][x][y];
}
int GraphMetricAntipode(int x)
{
    return x > n[type] ? x - n[type]: x + n[type];
}

void Load()
{
    printf("Loading metric\n");

    scanf("%d%d", &n[type], &m[type]);
    int N = n[type], M = m[type];
    for (int i = 1; i <= N; ++i)
        for (int j = 1; j <= N; ++j)
            dis[type][i][j] = (i == j? 0: 1e9);
    for (int i = 1; i <= M; ++i)
    {
        int u, v;
        double l;
        scanf("%d%d%lf", &u, &v, &l);
        dis[type][u][v] = dis[type][v][u] = min(1.0 * l, dis[type][u][v]);
    }
    for (int i = 1; i <= N; ++i)
        for (int j = 1; j <= N; ++j)
            for (int k = 1; k <= N; ++k)
                dis[type][j][k] = min(dis[type][j][k], dis[type][j][i] + dis[type][i][k]);
    
    double diameter = 0;
    for (int i = 1; i <= N; ++i)
        for (int j = 1; j <= N; ++j)
            diameter = max(diameter, dis[type][i][j]);
    printf("Diameter = %.3lf\n", diameter);
    for (int i = 1; i <= N; ++i)
        for (int j = 1; j <= N; ++j)
        {
            dis[type][i + N][j + N] = dis[type][i][j];
            dis[type][i + N][j] = dis[type][i][j + N] = diameter * 1.5 - dis[type][i][j];
        }
    printf("Construction finished\n");
}

int Edge[LENGTH];

WorkFunctionSimulator<3, int> GraphSimulator1(GraphMetricDistance), _GraphSimulator1(GraphMetricDistance), GraphSimulator2(GraphMetricDistance), _GraphSimulator2(GraphMetricDistance);

template<class T>
double verify(const T& WFA_Simulator1, const T& WFA_Simulator2) // for canonical potential
{
    // assert(K == 3);
    auto r = *WFA_Simulator2.requests.rbegin(), barr = WFA_Simulator2.antipode(r);
    auto config = Config(barr, barr, barr);
    double extended_cost = WFA_Simulator2.workFunction(config) - WFA_Simulator1.workFunction(config);
    // printf("extended cost = %.0lf\n", extended_cost);
    // printf("original potential = %.0lf\n", WFA_Simulator1.Potential_New(Edge));
    // printf("current potential = %.0lf\n", WFA_Simulator2.Potential_New(Edge));
    return WFA_Simulator2.Potential_New(Edge) - WFA_Simulator1.Potential_New(Edge) + eps >= extended_cost;
}

template<class T>
double verify(const T& WFA_Simulator) // for restricted canonical potential
{
    // assert(K == 3);
    auto r = *WFA_Simulator.requests.rbegin();
    return WFA_Simulator.Potential_New(r, Edge) - eps < WFA_Simulator.Potential_New(Edge);
}


void init()
{
    time_t seed = time(0);
    cout << "seed = " << seed << endl;
    srand(seed);
}

void LoadMetrics()
{
    type = 0;
    Load();
    n[0] /= 2;
    type = 1;
    Load();
}

void LoadExamples()
{
    CycleSimulator.antipode = CycleMetricAntipode;
    CycleSimulator.print = print<double>;
    CycleLength = 8;
    CycleSimulator.enumerator.clear();
    for (double i = 0; i < CycleLength; i += 1) CycleSimulator.enumerator.push_back(i);
    CycleSimulator.supp.clear();
    CycleSimulator.requests.clear();
    CycleSimulator.requests.push_back(5.);
    CycleSimulator.supp.push_back(make_pair(Config(5., 1., 6.), 9));
    CycleSimulator.supp.push_back(make_pair(Config(5., 1., 7.), 9));
    CycleSimulator.supp.push_back(make_pair(Config(5., 2., 3.), 10));
    CycleSimulator.supp.push_back(make_pair(Config(5., 2., 4.), 10));
    CycleSimulator.supp.push_back(make_pair(Config(5., 3., 6.), 8));
    CycleSimulator.supp.push_back(make_pair(Config(5., 3., 7.), 8));
    CycleSimulator.supp.push_back(make_pair(Config(5., 4., 6.), 8));
    CycleSimulator.supp.push_back(make_pair(Config(5., 4., 7.), 8));
    _CycleSimulator = CycleSimulator;
    _CycleSimulator.append(4.);
    
    type = 0;
    GraphSimulator1.antipode = GraphMetricAntipode;
    GraphSimulator1.print = print<int>;
    GraphSimulator1.enumerator.clear();
    for (int i = 1; i <= n[0] * 2; ++i) GraphSimulator1.enumerator.push_back(i);
    GraphSimulator1.supp.clear();
    GraphSimulator1.requests.clear();
    GraphSimulator1.requests.push_back(5);
    GraphSimulator1.supp.push_back(make_pair(Config(5, 2, 6), 9));
    GraphSimulator1.supp.push_back(make_pair(Config(5, 1, 2), 10));
    GraphSimulator1.supp.push_back(make_pair(Config(5, 3, 3), 11));
    GraphSimulator1.supp.push_back(make_pair(Config(5, 3, 4), 10));
    GraphSimulator1.supp.push_back(make_pair(Config(5, 3, 6), 9));
    GraphSimulator1.supp.push_back(make_pair(Config(5, 1, 3), 10));
    GraphSimulator1.supp.push_back(make_pair(Config(5, 4, 6), 8));
    GraphSimulator1.supp.push_back(make_pair(Config(5, 4, 1), 9));
    _GraphSimulator1 = GraphSimulator1;
    _GraphSimulator1.append(4);

    // Counterexample for Potential_3Circle
    type = 1;
    GraphSimulator2.antipode = GraphMetricAntipode;
    GraphSimulator2.print = print<int>;
    GraphSimulator2.enumerator.clear();
    for (int i = 1; i <= n[1] * 2; ++i) GraphSimulator2.enumerator.push_back(i);
    GraphSimulator2.supp.clear();
    GraphSimulator2.requests.clear();
    GraphSimulator2.requests.push_back(5);
    GraphSimulator2.supp.push_back(make_pair(Config(5, 6, 7), 0));
    GraphSimulator2.supp.push_back(make_pair(Config(5, 7, 8), 0));
    GraphSimulator2.supp.push_back(make_pair(Config(5, 8, 6), 0));
    GraphSimulator2.supp.push_back(make_pair(Config(5, 5, 6), 0));
    GraphSimulator2.supp.push_back(make_pair(Config(5, 5, 7), 0));
    GraphSimulator2.supp.push_back(make_pair(Config(5, 5, 8), 0));
    GraphSimulator2.supp.push_back(make_pair(Config(5, 5, 5), 0));
    _GraphSimulator2 = GraphSimulator2;
    _GraphSimulator2.append(1);
}

bool RandomRequest(int Round)
{
    CycleSimulatorR.antipode = CycleMetricAntipode;
    CycleSimulatorR.print = print<double>;
    CycleSimulatorR.enumerator.clear();
    for (double i = 0; i < CycleLength; i += 1) CycleSimulatorR.enumerator.push_back(i);
    CycleSimulatorR.init(Config(0., 0., 0.));
    CycleSimulatorR.requests.push_back(0.);
    double potential = CycleSimulatorR.Potential_New(Edge);
    for (int i = 1; i <= Round; ++i)
    {
        double r = rand() % 8, barr = CycleMetricAntipode(r);
        auto config = Config(barr, barr, barr);
        double wf = CycleSimulatorR.workFunction(config);
        CycleSimulatorR.append(r);
        double newPotential = CycleSimulatorR.Potential_New(Edge);
        if (newPotential - potential + eps < CycleSimulatorR.workFunction(config) - wf) return false;
        potential = newPotential;
    }
    return true;
}

int Start[LENGTH];
bool le(int A[], int B[])
{
    int p = 0;
    while (p < LENGTH && A[p] == B[p]) ++p;
    return p < LENGTH && A[p] < B[p];
}
bool leq(int A[], int B[])
{
    int p = 0;
    while (p < LENGTH && A[p] == B[p]) ++p;
    return p == LENGTH || A[p] < B[p];
}

vector<vector<int> > perms;
ofstream file("a.out");

void dfs(int x, int y)
{
    if (x == LENGTH - 1)
    {
        if (le(Edge, Start)) return;
        //if (le(Start, Edge)) return;
        static int _Edge[LENGTH];
        for (auto p: perms)
        {
            memset(_Edge, 0, sizeof(_Edge));
            for (int i = 0; i < LENGTH; ++i)
                for (int j = 0; j < LENGTH; ++j)
                    if ((Edge[i] >> j & 1))
                        _Edge[p[i]] |= 1 << p[j];
            if (le(_Edge, Edge)) return;
        }
        ofstream fout("config.txt");
        for (int i = 0; i < LENGTH; ++i) fout << bitset<LENGTH>(Edge[i]) << endl;
        fout << "Testing\n" << endl;
        fout.close();
        //LoadExamples();
        type = 0;
        if (!verify(_GraphSimulator1)) return;
        //if (!verify(GraphSimulator1, _GraphSimulator1)) return;
        type = 1;
        if (!verify(_GraphSimulator2)) return;
        //if (!verify(GraphSimulator2, _GraphSimulator2)) return;
        if (!verify(_CycleSimulator)) return;
        //if (!verify(CycleSimulator, _CycleSimulator)) return;
        //if (!RandomRequest(100)) return;
        for (int i = 0; i < LENGTH; ++i) file << bitset<LENGTH>(Edge[i]) << endl;
        file << "Passed All Tests\n" << endl;
        for (int i = 0; i < LENGTH; ++i) cout << bitset<LENGTH>(Edge[i]) << endl;
        cout << "Passed All Tests\n" << endl;
        return;
    }
    int nx = x, ny = y - 1;
    if (y == x + 1) nx = x + 1, ny = LENGTH - 1;
    dfs(nx, ny);
    Edge[x] ^= 1 << y; Edge[y] ^= 1 << x;
    dfs(nx, ny);
    Edge[x] ^= 1 << y; Edge[y] ^= 1 << x;
}

void InitPermutations()
{
    vector<int> e;
    for (int i = 0; i < LENGTH; ++i) e.push_back(i);
    vector<int> p3;
    for (int i = 0; i < 3; ++i) p3.push_back(i);
    for (int i = 0; i < 2; ++i)
    {
        for (int i = 0; i < 2; ++i)
        {
            for (int i = 0; i < 2; ++i)
            {
                do{
                    vector<int> p;
                    for (int i = 0; i < 3; ++i)
                    {
                        p.push_back(e[p3[i] * 2]);
                        p.push_back(e[p3[i] * 2 + 1]);
                    }
                    p.push_back(6);
                    perms.push_back(p);
                }while (next_permutation(p3.begin(), p3.end()));
                swap(e[4], e[5]);
            }
            swap(e[2], e[3]);
        }
        swap(e[0], e[1]);
    }
}

int main()
{
    init();
    LoadMetrics();
    InitPermutations();
    LoadExamples();

    ifstream fin("config.txt");
    for (int i = 0; i < LENGTH; ++i)
    {
        bitset<LENGTH> B;
        fin >> B;
        Start[i] = B.to_ulong();
    }

    dfs(0, LENGTH - 1);

    fin.close();
    return 0;
}
