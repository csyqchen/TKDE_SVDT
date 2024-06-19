#include <iostream>
#include <queue>
#include <vector>
#include <unordered_set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <set>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <iomanip>

using namespace std;

vector <int>* t;
int tree_n;
double* feq1, * feq2 , mu;
int* parent;
set <int> ans, ans1, ans2;
vector <int> S, S1, S2;
int K, K1, beta;
double *marginSim, *marginDif;
vector<vector<int>> dist;
vector<vector<int>> desc;

int* RoundSim, * RoundDif;

struct Candidate {
    int id;
    double similarity;
    double difference;
    double metricSim, metricDif;

    bool operator<(const Candidate& other) const {
        return std::max(similarity, difference) < std::max(other.similarity, other.difference);
    }
};

using PQ = std::priority_queue<Candidate>;

PQ candidateByScore;

struct DistriSim{
    int IDsim;
    double Vsim;
};

struct DistriDif{
    int IDdif;
    double Vdif;
};

struct CompareDistriSim{
    bool operator()(const DistriSim& d1, const DistriSim& d2){
        return d1.Vsim < d2.Vsim;
    }
};

struct CompareDistriDif{
    bool operator()(const DistriDif& d1, const DistriDif &d2){
        return d1.Vdif < d2.Vdif;
    }
};

DistriSim dictSim;
DistriDif dictDif;

std::vector<std::priority_queue<DistriSim, std::vector<DistriSim>, CompareDistriSim>> DistriSimQueue(tree_n);
std::vector<std::priority_queue<DistriDif, std::vector<DistriDif>, CompareDistriDif>> DistriDifQueue(tree_n);

vector<double> Normlization(std::vector<double> array, bool y) {
    std::vector<double> result;

    double sum = 0.0;
    for (int i = 0; i < array.size(); i++) {
        sum += array[i];
    }
    
    if (sum < 1e-6) {
        for (int i = 0; i < array.size(); i++)
            result.push_back(1.0 / array.size());
        }
    else {
        for (int i = 0; i < array.size(); i++)
            result.push_back(array[i] / sum);
        }
    return result;
}
    
struct MinValue {
    double value;
    int index;
};

bool compareMinValue(const MinValue& a, const MinValue& b) {
    return a.value > b.value;
}

double get_mu() {
    vector<MinValue> minValues;
    for (int i = 1; i <= tree_n; i++) {
        MinValue minValue;
        minValue.value = min(feq1[i], feq2[i]);
        minValue.index = i;
        minValues.push_back(minValue);
    }

    sort(minValues.begin(), minValues.end(), compareMinValue);

    double sumsim = 0.0;
    double sumdif = 0.0;
    for (int i = 0; i < tree_n; i++) {
        sumsim += minValues[i].value;
        sumdif += abs(feq1[minValues[i].index] - feq2[minValues[i].index]);
    }

    return sumsim / sumdif;
}

double hellingerDistance(int x)
{
    vector<double> p, q;

    while ((!DistriSimQueue[x].empty()) and (!DistriDifQueue[x].empty())){
        int pvalue1 = feq1[DistriSimQueue[x].top().IDsim];
        int pvalue2 = feq2[DistriSimQueue[x].top().IDsim];
        int qvalue1 = feq1[DistriDifQueue[x].top().IDdif];
        int qvalue2 = feq2[DistriDifQueue[x].top().IDdif];
        p.push_back(pvalue1); p.push_back(qvalue1);
        q.push_back(pvalue2); q.push_back(qvalue2);
        DistriSimQueue[x].pop();
        DistriDifQueue[x].pop();
    }
    
    p = Normlization(p, true);
    q = Normlization(q, false);
    double distance = 0.0;
    for (int i = 0; i < min(p.size(), q.size()); i++) {
        double sqrt_p = sqrt(p[i]);
        double sqrt_q = sqrt(q[i]);
        distance += (sqrt_p - sqrt_q) * (sqrt_p - sqrt_q);
    }
    distance = sqrt(distance) / sqrt(2.0);
    return distance;
}

double Summary_Impact(double frequency1, double frequency2, int level, bool IsSim, int x)
{
    double temp, score;
    if ((frequency1 >= 1e-6) || (frequency2 >= 1e-6)) { temp = min(frequency1, frequency2) / max(frequency1, frequency2); }
    else { temp = 0; }
    if (IsSim) {
        score = max(frequency1, frequency2) * temp / sqrt(level);
    } else {
        score = max(frequency1, frequency2) * mu *(1 - temp) / sqrt(level);
    }
    return score;
}


double Margin_DFS(int x, int level, bool IsSim)
{
    double sum = 0;
    set<int> visited;
    if (IsSim ? ans1.find(x) == ans1.end():ans2.find(x) == ans2.end())
    {
        visited.insert(x);
        sum += Summary_Impact(feq1[x], feq2[x], level, IsSim, x);
    }
    level++;
    for (int i = 0; i < t[x].size(); i++)
    {
        int y = t[x][i];
        if ((IsSim ? ans1.find(x) == ans1.end():ans2.find(x) == ans2.end()) && visited.find(y) == visited.end())
        {
            visited.insert(y);
            sum += Margin_DFS(y, level, IsSim);
        }
    }
    return sum;
}

double Margin_Update(int x, int level, bool IsSim)
{
    double sum = Margin_DFS(x, level, IsSim);
    double reduce = 0;
    
    int p = x;
    int diff_level = 0;
    while (p != -1 && (IsSim ? ans1.find(p) == ans1.end():ans2.find(p) == ans2.end()))
    {
        p = parent[p];
        diff_level++;
    }
    if (p != -1)
    {
        reduce = Margin_DFS(x, level + diff_level, IsSim);
    }
    return sum - reduce;
}

void Margin_Init(int n)
{
    Candidate temp;
    int times = 0;
    beta = 50;
    DistriSimQueue.resize(tree_n);
    DistriDifQueue.resize(tree_n);
    dist.resize(n + 1);
    desc.resize(n + 1);
    for (int i = 1; i <= n; i++)
    {
        times = 1;
        if ((feq1[i] >= 1e-6) || (feq2[i] >= 1e-6))
        {
            int p = i;
            int level = 1;
            while (p != -1)
            {
                marginSim[p] += Summary_Impact(feq1[i], feq2[i], level, true, p);
                marginDif[p] += Summary_Impact(feq1[i], feq2[i], level, false, p);
                if (!t[p].empty() ){
                    if (DistriSimQueue[p].size() < beta) {
                        DistriSim temp;
                        temp.IDsim = i;
                        temp.Vsim = min(feq1[i], feq2[i]);
                        DistriSimQueue[p].push(temp);
                    }
                    else{
                        if (min(feq1[i], feq2[i]) <= DistriSimQueue[p].top().Vsim)
                        {}
                        else {
                            DistriSim temp = DistriSimQueue[p].top();
                            DistriSimQueue[p].pop();
                            temp.IDsim = i;
                            temp.Vsim = min(feq1[i], feq2[i]);
                            DistriSimQueue[p].push(temp);
                        }
                    }
                    if (DistriDifQueue[p].size() < beta){
                        DistriDif temp1;
                        temp1.IDdif = i;
                        temp1.Vdif = abs(feq1[i] - feq2[i]);
                        DistriDifQueue[p].push(temp1);
                    }
                    else{
                        if (abs(feq1[p] - feq2[p]) <= DistriDifQueue[p].top().Vdif)
                        {}
                        else {
                            DistriDif temp1 = DistriDifQueue[p].top();
                            DistriDifQueue[p].pop();
                            temp1.IDdif = i;
                            temp1.Vdif = abs(feq1[i] - feq2[i]);
                            DistriDifQueue[p].push(temp1);
                        }
                    }
                }
                if (parent[p] != -1){
                    dist[parent[p]].push_back(level);
                    desc[parent[p]].push_back(p);
                }
                level++;
                p = parent[p];
                times++;
            }
        }
    }

    for (int i = 1; i <= n; i++)
    {
        temp.id = i;
        if (t[i].size()!= 0 ) {
            double temp2 = hellingerDistance(i);
            temp.similarity = marginSim[i] * (1 - temp2);
            if (temp2 != 0){
                temp.difference = marginDif[i] * temp2;
            }
            else {
                temp.difference = marginDif[i];
            }
        }
        else {
            temp.similarity = marginSim[i];
            temp.difference = marginDif[i];
        }
        candidateByScore.push(temp);
        RoundSim[i] = 0;
        RoundDif[i] = 0;
        printf("#%d Similarity:%.2lf \n", i, marginSim[i]);
        printf("#%d Difference:%.2lf \n", i, marginDif[i]);
    }
}
    
double GSDVT(int K){
        int k = K;
        bool isSim = true;
        double impact = 0.0, diversity = 0.0;
        ans.clear();
        ans1.clear();
        ans2.clear();
        S.clear();
        S1.clear();
        S2.clear();

        while (S1.size() + S2.size() < k && !candidateByScore.empty()) {
            
            Candidate x = candidateByScore.top();
            candidateByScore.pop();
            
            if (x.similarity >= x.difference) {
                isSim = true;
                if (RoundSim[x.id] < S1.size() && S1.size() != 0 ) {
                    x.similarity = Margin_Update(x.id, 1, true);
                    isSim = true;
                    candidateByScore.push(x);
                    RoundSim[x.id] = S1.size();
                    printf("Sim Margin: %lf \n", x.similarity);
                    continue;
                }
            } else {
                isSim = false;
                if (RoundDif[x.id] < S2.size() && S2.size() != 0 ) {
                    x.difference = Margin_Update(x.id, 1, false);
                    isSim = false;
                    candidateByScore.push(x);
                    RoundDif[x.id] = S2.size();
                    printf("Dif Margin: %lf \n", x.difference);
                    continue;
                }
            }
            
            S.push_back(x.id);
            isSim ? S1.push_back(x.id):S2.push_back(x.id);
            ans.insert(x.id);
            isSim ? ans1.insert(x.id):ans2.insert(x.id);
        }
        
        printf("S1: %d S2: %d\n", S1.size(), S2.size());
        return 0;
}

int main()
{
    FILE* fin, * fout;
    int n;
    int x, y, z;

    fin = fopen("./graph.txt", "r");

    int m;

    fscanf(fin, "%d %d", &m, &n);
    printf("%d %d\n", m, n);
    tree_n = n;
    beta = 50;
    t = new vector <int>[n + 1];
    feq1 = new double[n + 1];
    feq2 = new double[n + 1];
    marginSim = new double[n + 1];
    marginDif = new double[n + 1];
    parent = new int[n + 1];
    RoundSim = new int[n + 1];
    RoundDif = new int[n + 1];
    cover = new bool[n + 1];
    metricSim = new double [n + 1];
    metricDif = new double [n + 1];

    ans.clear();
    memset(feq1, 0, sizeof(double) * (n + 1));
    memset(feq2, 0, sizeof(double) * (n + 1));
    memset(parent, -1, sizeof(int) * (n + 1));
    memset(marginSim, 0, sizeof(double) * (n + 1));
    memset(marginDif, 0, sizeof(double) * (n + 1));
    memset(RoundSim, 0, sizeof(int) * (n + 1));
    memset(RoundDif, 0, sizeof(int) * (n + 1));
    memset(cover, false, sizeof(bool) * (n + 1));
    memset(metricSim, 0, sizeof(double) * (n + 1));
    memset(metricDif, 0, sizeof(double) * (n + 1));

    
    for (int i = 0; i < m - 1 ; i++)
    {
        fscanf(fin, "%d %d", &x, &y);
        
        if (parent[y] == -1)
        {
            t[x].push_back(y);
            parent[y] = x;
        }
    }
    fclose(fin);

    fin = fopen("./freq1.txt", "r");

    fscanf(fin, "%d", &n);
    printf("%d\n", n);
        
    for (int i = 0; i < n; i++)
    {
        fscanf(fin, "%d %d", &x, &y);
        feq1[x] = y;
    }
    fclose(fin);
    
    fin = fopen("./freq2.txt", "r");
 
    fscanf(fin, "%d", &n);

    for (int i = 0; i < n; i++)
    {
        fscanf(fin, "%d %d", &x, &y);
        feq2[x] = y;
    }
    fclose(fin);
    
    K = 10;
    
    mu = get_mu();
    
    Margin_Init(tree_n);
    GSDVT(K);
    
    for (int k=90; k<=90; k+=10)
    {
        ans.clear();
        for (int i = 0; i < k && i < S.size(); i++)
        {
            ans.insert(S[i]);
        }
        for (int i = 0; i < k && i < S1.size(); i++)
        {
            ans1.insert(S1[i]);
        }
        for (int i = 0; i < k && i < S2.size(); i++)
        {
            ans2.insert(S2[i]);
        }
        printf("#%d\n", (int)ans.size());
    }

    for (int i = 0; i < K; i++)
    {
        printf("#%d ", S[i]);

        int x = S[i];

        while (1)
        {
            printf("%d", x);
            x = parent[x];
            if (x == -1) break;
            printf("-->");
        }
        printf("\n");
    }
    
    for (int i = 0; i < ans1.size(); i++)
    {
        printf("#%d Sim:", S1[i]);

        int x = S1[i];

        while (1)
        {
            printf("%d", x);
            x = parent[x];
            if (x == -1) break;
            printf("-->");
        }
        printf("\n");
    }
    
    for (int i = 0; i < ans2.size(); i++)
    {
        printf("#%d Dif:", S2[i]);

        int x = S2[i];

        while (1)
        {
            printf("%d", x);
            x = parent[x];
            if (x == -1) break;
            printf("-->");
        }
        printf("\n");
    }
    printf("%lf\n", mu);
    return 0;
}
