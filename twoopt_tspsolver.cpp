#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>

#define MODE 2	//1 : two-opt
				//2 : greedy -> tow-opt
#define RANK 10
#define EMPTY -1

using namespace std;

typedef struct{
	double x;
	double y;
}Point;

double Dis(Point a, Point b)
{
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	return sqrt(dx*dx + dy*dy);
}

typedef struct{
	int num;
	double length;
}Distance;

bool operator<(const Distance& left, const Distance& right)
{
	return left.length < right.length;
}


void Flip(const int b, const int d, vector< int > *tour, const vector< Point > &points)
{

	vector< int > top, middle, under;
	vector< int >::iterator it = (*tour).begin();


	for (; *it != b; ++it)
		top.push_back(*it);

	for (; *it != d; ++it){
		middle.push_back(*it);

		//dがbよりも手前に存在する場合
		//[head ... c][d ... a][b ... end] -> [d ... a][c ... head end ... b]

		if (it+1 == (*tour).end()){

			it = (*tour).begin();
			for (; *it != d; ++it)
				middle.push_back(*it);

			top.clear();
			for (; *it != b; ++it)
				top.push_back(*it);

			reverse(middle.begin(), middle.end());
			top.insert(top.end(), middle.begin(), middle.end());
			(*tour) = top;
			return;
		}
	}
	for (; it != (*tour).end(); ++it)
		under.push_back(*it);


	reverse(middle.begin(), middle.end());

	top.insert(top.end(), middle.begin(), middle.end());
	top.insert(top.end(), under.begin(), under.end());
	(*tour) = top;
}

inline bool isAllTrue(vector<bool> &vec)
{
	for (auto v : vec){
		if (v == false)return false;
	}
	return true;
}

void greedy(vector< int > *tour, double *tlength, vector< Point > &points, int n)
{
	vector< bool > chkTour;
	for (int i = 0; i < n; i++){
		chkTour.push_back(false);
	}

	int now = 0;
	(*tour).push_back(now);
	chkTour[now] = true;

	while(!isAllTrue(chkTour)){

		//nowから各点までの距離と番号のデータ(Distance)をpointsLengthに格納
		vector< Distance > pointsLength;
		for (int i = 0; i < n; i++){
			if (chkTour[i] == false){
				Distance tmp;
				tmp.length = Dis(points[now], points[i]);
				tmp.num = i;
				pointsLength.push_back(tmp);

			}
		}
		//lengthでソート
		sort(pointsLength.begin(), pointsLength.end());

		/*
		//for debug
		for (auto p : pointsLength){
			cout << p.num << " ";
		}
		cout << endl;
		*/

		//nowから距離が最短の点を選び、次のnowとする。
		now = pointsLength[0].num;
		(*tour).push_back(now);
		*tlength += pointsLength[0].length;
		chkTour[now] = true;
	}

}

double solver(vector< Point > &points, int n)
{
	//初期巡回路
	vector< int > tour;
	double length = 0;
	//

	switch (MODE){
		case 1 : for (int i = 0; i < n; i++){
			tour.push_back(i);
			length += Dis(points[i], points[(i + 1) % n]);
		}
				 break;
		case 2 : 
		default: greedy(&tour, &length, points, n);
			break;
	}
	
	
	two_opt_start:
	for (int i = 0; i < n-1; i++){

		int a = tour[i];
		int b = tour[(i + 1)%n];

		//bから最短の点をRANK個引っ張る
		int nearPoints[RANK];
		for (int j = 0; j < RANK; j++)nearPoints[j] = EMPTY;

		vector< Distance > pointsLength;

		//pointsLengthの取得
		Point pa = points[a];
		for (int i = 0; i < n; i++){
			Distance tmp;
			tmp.num = i;
			tmp.length = Dis(pa, points[i]);
			if (tmp.length == 0)continue;
			pointsLength.push_back(tmp);
		}

		sort(pointsLength.begin(), pointsLength.end());

		for (int j = 0; j < RANK; j++){
			if (j > n)continue;
			nearPoints[j] = pointsLength[j].num;
		}

		pointsLength.clear();
		
		/*
		//for debug
		for (auto np : nearPoints){
			cout << np << " ";
		}
		cout << endl;
		*/

		for (auto np : nearPoints){
			//cの確定
			int c = np;
			if (a == c || b == c || c == EMPTY || c == n-1)continue;
			vector< int >::iterator it = tour.begin();

			//dの確定
			while (*it != c)it++;
			int d;
			if (it+1 == tour.end())d = *tour.begin();
			else d = *(it + 1);
			if (a == d)continue;
			Point pb = points[b], pc = points[c], pd = points[d];

			double tmp = Dis(pa, pb) + Dis(pc, pd) - Dis(pa, pc) - Dis(pb, pd);

			if (tmp > 0.01){

				//cout << "tmp = " << tmp << endl;
				//cout << a << "," << b << "," << c << "," << d << endl;

				Flip(b, d, &tour, points);
				length -= tmp;

				
				//for debug
				//for (auto t : tour){cout << t << "->";}
				//cout << endl;
				

				cout << "length = " << length << endl;
				goto two_opt_start;
			}	
		}
	}
	for (auto t : tour){cout << t << "->";}
	cout << endl;
	return length;
}

int main()
{
	int n;
	cin >> n;
	
	vector<Point> points;

	for(int i = 0; i < n; i++){
		Point tmp;
		cin >> tmp.x >> tmp.y;
		//cout << tmp.x << tmp.y << endl;

		points.push_back(tmp);
	}

	time_t start, end;

	start = clock();
	cout << "length = " << solver(points, n) << endl;
	end = clock();

	cout << end - start << "sec" << endl;
	return 0;
}