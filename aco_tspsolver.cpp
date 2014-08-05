#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iomanip>

//Define���ꂽ�萔
#define MODE 2			//1 : AS0
						//2 : ASrank
						//3 : MIN-MAX
#define TAU_DEFAULT 0.1	//����pheromon�̗�
#define ANTS_RATE 0.5	//�_�̐��ɑ΂��Ĕz�u�����A���̔䗦
#define ALPHA 1			//pheromon�̉e���x
#define BETA 5			//�����̉e���x
#define RHO 0.5			//pheromon�̏�����
#define Q 100.0			//��������pheromon�̔{��
#define LOOP_MAX 100	//�ő��LOOP��
#define MIN_PHERO 0.001	//pheromon�̍ŏ��l

#define ELITE_RATE 10	//elite�A���̐�

//�O���[�o���ϐ�(num�����͎���Ƃ���ŗp����̂�)
int NUM;

//����
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

//�O�p�^�z����ꎟ���z��Ƃ��Ĉ������߂̎�
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

//�������v�Z
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
	Point *points;							//�_�̏��
//	Node *nodes = new Node[getTrinum(NUM, NUM)];
	double *phero = new double[getTrinum(NUM, NUM)];

											//2�_�Ԃ̏��
	std::vector<Ant> ants;					//�A����vector
	std::vector<Ant> elite;					//�D�G�ȃA��
	int num_elite;

public:
	Towns();
	Towns(Point *p);
	void setAnts();
	double getDis(int i, int j) const;		//���������߂�
	double getPhero(int i, int j) const;	//pheromon�����߂�
	int getWhere(Ant &a);					//�������֍s�������߂�
	Towns& gotoPoint(Ant &a, int next);		//Ant a��_now����_next�ֈړ�������
	void moveAnts();						//�A�����ړ�
	Towns& setPhero(int i, int j, double p);
											//pheromon��ݒ�
	void renewPhero();						//pheromon���X�V

	Towns& setRank();						//elite�A���̌���
	bool isSolved(){
		if (MODE == 2)
			return int(elite[0].length) == int(elite[num_elite - 1].length);
		else return false;
	};
									//�g�b�v�̃G���[�g�A���ƍŉ��ʂ̃G���[�g�A���������o�͂Ȃ�Β�~

	double getBestLength(){ return elite[0].length; };
	double getAverageLength();
	void Towns::showPhero();			//for debug
	void Towns::showAntsRoot();			//for debug
	~Towns(){ delete[] phero; }
};

Towns::Towns(Point *p){
	points = p;
	num_elite = ELITE_RATE;

	//nodes��ݒ�
	for (int i = 0; i < NUM; ++i){
		for (int j = 0; j <= i; ++j){
			phero[getTrinum(i, j)] = TAU_DEFAULT;
		}
	}
}

void Towns::setAnts()
{
	//ants��������
	ants.clear();

	//ants��ݒ�
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

//�܂��K��Ă��Ȃ��_�̏��
typedef struct{
	int point;
	double pij;
}notVisited;

int Towns::getWhere(Ant &a)
{
	//i = a.now
	//p(i, j) = pij/psum = pow(phero(i, j), ALPHA) * pow(1.0/dis(i, j), BETA)
	//			Sigma(ForAll k [visited[k] == true])(pow(phero(i, k), ALPHA) * pow(1/dis(i, k), BETA))


	//����start�n�_�ȊO�̓_��S�ĉ�����Ȃ�΁Astart�n�_�ɋA��
	bool bl = true;
	for (int j = 0; j < NUM; ++j){
		bl = bl && a.visited[j];
	}
	if (bl)return a.root[0];

	std::vector<notVisited> nVisited;
	double psum = 0;
	int i = a.now;
	for (int j = 0; j < NUM; ++j){

		//�܂��K��Ă��Ȃ��ꏊ�݂̂�I��
		if (a.visited[j] == false){

			notVisited nv;
			double tmp = pow(getPhero(i, j), ALPHA) * pow(1.0 / getDis(i, j), BETA);

			psum += tmp;
			nv.point = j;
			nv.pij = tmp;
			nVisited.push_back(nv);
		}
	}

	//pij���m�����z�̔z��Ƃ���
	double probSum = 0;	//�m���̘a
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

	//sigmaTau�̌v�Z
	switch (MODE){
		//AS
	case 1: {
				//elite�̐ݒ�
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

				//elite�̐ݒ�
				elite.clear();
				setRank();

				//elite�̒ʉ߂������̂�pheromon���ǉ������
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

	//phero�̍X�V
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

//�ȉ�debug�p
//���݂�pheromon�ʂ�\��
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

//�A���̋O�Ղ�\��
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

	//���͌`����
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
		//�񐔋L�^
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
