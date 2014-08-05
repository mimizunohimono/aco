#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iomanip>

//Defineされた定数
#define MODE 2			//1 : AS0
						//2 : ASrank
						//3 : MIN-MAX
#define TAU_DEFAULT 0.1	//初期pheromonの量
#define ANTS_RATE 0.5	//点の数に対して配置されるアリの比率
#define ALPHA 1			//pheromonの影響度
#define BETA 5			//距離の影響度
#define RHO 0.5			//pheromonの蒸発率
#define Q 100.0			//増加するpheromonの倍率
#define LOOP_MAX 100	//最大のLOOP回数
#define MIN_PHERO 0.001	//pheromonの最小値

#define ELITE_RATE 10	//eliteアリの数

//グローバル変数(numだけは至るところで用いるので)
int NUM;

//乱数
inline void InitRand()
{
	srand((unsigned int)time(NULL));
}

inline bool getRand(double rate)
{
	if (rand() / (RAND_MAX + 1.0) < rate)
		return true;
	else
		return false;
}

//三角型配列を一次元配列として扱うための式
inline int getTrinum(int i, int j)
{
	if (i > j)
		return i * (i + 1) / 2 + j;
	else return j * (j + 1) / 2 + i;
}

typedef struct{
	double x;
	double y;
}Point;

//距離を計算
double Dis(Point a, Point b)
{
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	return double(sqrt(dx*dx + dy*dy));
}

typedef struct{
	double dis;
	double phero;
}Node;

typedef struct{
	std::vector<int> root;
	bool *visited;
	double length;
	int now;
}Ant;

class Towns
{
	Point *points;							//点の情報
//	Node *nodes = new Node[getTrinum(NUM, NUM)];
	double *phero = new double[getTrinum(NUM, NUM)];

											//2点間の情報
	std::vector<Ant> ants;					//アリのvector
	std::vector<Ant> elite;					//優秀なアリ
	int num_elite;

public:
	Towns();
	Towns(Point *p);
	void setAnts();
	double getDis(int i, int j) const;		//距離を求める
	double getPhero(int i, int j) const;	//pheromonを求める
	int getWhere(Ant &a);					//次何処へ行くか求める
	Towns& gotoPoint(Ant &a, int next);		//Ant aを点nowから点nextへ移動させる
	void moveAnts();						//アリを移動
	Towns& setPhero(int i, int j, double p);
											//pheromonを設定
	void renewPhero();						//pheromonを更新

	Towns& setRank();						//eliteアリの決定
	bool isSolved(){
		if (MODE == 2)
			return int(elite[0].length) == int(elite[num_elite - 1].length);
		else return false;
	};
									//トップのエリートアリと最下位のエリートアリが同じ出力ならば停止

	double getBestLength(){ return elite[0].length; };
	double getAverageLength();
	void Towns::showPhero();			//for debug
	void Towns::showAntsRoot();			//for debug
	~Towns(){ delete[] phero; }
};

Towns::Towns(Point *p){
	points = p;
	num_elite = ELITE_RATE;

	//nodesを設定
	for (int i = 0; i < NUM; ++i){
		for (int j = 0; j <= i; ++j){
			phero[getTrinum(i, j)] = TAU_DEFAULT;
		}
	}
}

void Towns::setAnts()
{
	//antsを初期化
	ants.clear();

	//antsを設定
	for (int i = 0; i < NUM; ++i){
		if (getRand(ANTS_RATE)){
			Ant tmp;
			tmp.root.push_back(i);
			tmp.visited = new bool[NUM];
			for (int j = 0; j < NUM; ++j)
				tmp.visited[j] = (i == j);
			tmp.length = 0;
			tmp.now = i;
			ants.push_back(tmp);
		}
	}
}

double Towns::getDis(int i, int j) const
{
	return Dis(points[i], points[j]);
	//return nodes[getTrinum(i, j)].dis;
}

double Towns::getPhero(int i, int j) const
{
	if (i == j)return 0;

	return phero[getTrinum(i, j)];
}

//まだ訪れていない点の情報
typedef struct{
	int point;
	double pij;
}notVisited;

int Towns::getWhere(Ant &a)
{
	//i = a.now
	//p(i, j) = pij/psum = pow(phero(i, j), ALPHA) * pow(1.0/dis(i, j), BETA)
	//			Sigma(ForAll k [visited[k] == true])(pow(phero(i, k), ALPHA) * pow(1/dis(i, k), BETA))


	//もしstart地点以外の点を全て回ったならば、start地点に帰る
	bool bl = true;
	for (int j = 0; j < NUM; ++j){
		bl = bl && a.visited[j];
	}
	if (bl)return a.root[0];

	std::vector<notVisited> nVisited;
	double psum = 0;
	int i = a.now;
	for (int j = 0; j < NUM; ++j){

		//まだ訪れていない場所のみを選択
		if (a.visited[j] == false){

			notVisited nv;
			double tmp = pow(getPhero(i, j), ALPHA) * pow(1.0 / getDis(i, j), BETA);

			psum += tmp;
			nv.point = j;
			nv.pij = tmp;
			nVisited.push_back(nv);
		}
	}

	//pijを確率分布の配列とする
	double probSum = 0;	//確率の和
	int count = 0;

	for (auto nv : nVisited){
		double tmp = nv.pij;
		tmp /= psum;
		probSum += tmp;
		nVisited[count++].pij = probSum;
	}

	double r = rand() / (RAND_MAX + 1.0);

	for (auto nv : nVisited){
		if (r <= nv.pij){
			return nv.point;
		}
	}
	return 0;
}

Towns& Towns::gotoPoint(Ant &a, int next)
{
	a.root.push_back(next);
	a.visited[next] = true;
	a.length += getDis(a.now, next);
	a.now = next;

	return *this;
}

void Towns::moveAnts()
{
	std::vector<Ant>::iterator it = ants.begin();

	while (it != ants.end()){
		int tmp = getWhere(*it);
		gotoPoint(*it, tmp);
		*it++;

	}
}

Towns& Towns::setPhero(int i, int j, double p)
{
	double tmp = p;
	if (i == j)tmp = 0;

	phero[getTrinum(i, j)] = tmp;
	return *this;
}

void Towns::renewPhero()
{
	//phero = (1 - RHO)*phero + sigma(dt(i, j))
	//dt(i, j) = (f(X^k).contain((i, j))) ? 1/f(X^k) 0

	std::vector<double> sigmaTau(getTrinum(NUM, NUM));

	for (int i = 0; i < getTrinum(NUM, NUM); ++i)sigmaTau[i] = 0;

	//sigmaTauの計算
	switch (MODE){
		//AS
	case 1: {
				//eliteの設定
				elite.clear();
				setRank();

				for (auto a : ants){
					int pre = a.root[0];
					for (auto r : a.root){
						sigmaTau[getTrinum(r, pre)] += Q / a.length;
						pre = r;
					}
				}
	}break;
		//ASrank
	case 2: {

				//eliteの設定
				elite.clear();
				setRank();

				//eliteの通過した道のみpheromonが追加される
				for (int i = 0; i < num_elite; ++i){
					int pre = elite[i].root[0];
					for (auto next : elite[i].root){

						//sigmaTau = Sigma(1, sig){(sig - mu) * Q / Lm}
						sigmaTau[getTrinum(next, pre)] +=
							double((num_elite - i) * Q) / elite[i].length;
						pre = next;
					}
				}
	}break;
		//MMAS
	case 3: {
				elite.clear();
				setRank();

				int pre = elite[0].root[0];
				
				for (auto next : elite[0].root){

					sigmaTau[getTrinum(next, pre)] +=
						Q / elite[0].length;
					pre = next;
				}
	}
	default: break;
	}

	//pheroの更新
	for (int i = 0; i < NUM; ++i){
		for (int j = 0; j < i; ++j){

	 		double tmp = (1 - RHO) * getPhero(i, j);
			tmp += sigmaTau[getTrinum(i, j)];

			if (MODE == 3 && tmp > Q / (elite[0].length * RHO))
				tmp = Q / (elite[0].length * RHO);
			if (tmp < MIN_PHERO)tmp = MIN_PHERO;
			setPhero(i, j, tmp);

		}
	}
}

bool comp(const Ant& left, const Ant& right) {
	return left.length < right.length;
}

Towns& Towns::setRank()
{
	std::vector<Ant> cpAnts = ants;
	sort(cpAnts.begin(), cpAnts.end(), comp);
	for (int i = 0; i < num_elite; ++i){
		elite.push_back(cpAnts[i]);
	}
	return *this;
}

double Towns::getAverageLength(){
	double sum = 0;
	for (auto e : elite){
		sum += e.length;
	}
	return sum / num_elite;
}

//以下debug用
//現在のpheromon量を表示
void Towns::showPhero()
{
	for (int i = 0; i < NUM; i++){
		for (int j = 0; j < NUM; j++){
			std::cout << std::fixed << std::setprecision(3) << getPhero(i, j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//アリの軌跡を表示
void Towns::showAntsRoot()
{
	for (auto a : ants){
		for (auto r : a.root)
			std::cout << r << " ";
		std::cout << " length = " << a.length << std::endl;
		std::cout << std::endl;
	}
}

int main(int argc, int argv[]){

	//入力形式は
	//n
	//x1 y1
	//x2 y2
	//...
	//xn yn

	InitRand();
	std::cin >> NUM;
	Point *points = new Point[NUM];
	double tmp1, tmp2, n;
	for (int i = 0; i < NUM; i++){
		std::cin >> n >> tmp1 >> tmp2;
		points[i].x = tmp1;
		points[i].y = tmp2;
	}

	//Step1 Initialize
	Towns *towns;
	towns = new Towns(points);

	bool tmpsolve;

	clock_t start, end;
	//towns->showPhero();
	int k = 0;
	start = clock();
	do{
		//回数記録
		//std::cout <<"k = " << k << std::endl;
		towns->setAnts();

		//Step2 Move ants
		for (int i = 0; i < NUM; i++){
			towns->moveAnts();
		}

		//Step3 Renew phero
		towns->renewPhero();
		//towns->showPhero();

		tmpsolve = towns->isSolved();
		
	} while (!tmpsolve && k++ < LOOP_MAX);

	end = clock();
	//towns->showPhero();
//	towns->showAntsRoot();
	
	//Output
	std::cout << "best solution == " << towns->getBestLength() << std::endl;
	std::cout << "average solution == " << towns->getAverageLength() << std::endl;
	std::cout << (double)(end - start) / CLOCKS_PER_SEC << "sec" << std::endl;
	//Finish
	delete points;
	delete towns;
	return 0;
}
