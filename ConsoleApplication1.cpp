// ConsoleApplication1.cpp : �������̨Ӧ�ó������ڵ㡣
//
#include "stdafx.h"
#include<math.h>
#include<iostream>
#include <cstring>
#include<time.h>
#include<string.h>
#include<stdio.h>
#include<list>
#include <stdlib.h>
#include <algorithm>
#include <cstdio>  
#include <cstdlib>  
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include<conio.h> 
#include<time.h>
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

using namespace std;
#define M_A_X 9999999999999999999  //MAX����Ϊ30��
typedef vector< vector<int> >  vvt;
#ifndef Index
struct Index
{
	int x, y;
};
#endif

vector<Index> test(int inN ,int ink ,double inr )
{

	vector<Index> mytest;
	double r = inr;
	int TIMES = 1;
	int N = inN;
	int k = ink;
	while (TIMES--)
	{
		printf("%d\n", N);

		vvt  vv, vv1, matrix, bit, bit1;
		vvt::iterator iter;
		vector<int>::iterator intiter, intiter2, intiter3, intiter7;
		vector<double>::iterator intiter8;
		vector<int> vi, vi1, row, row1, row2, vi2, vi3, s, breaknode, ND, ND1, ND2, nextposition, ND4, ND5;
		vector<double> capacity, ND3, ND6, ND7, B;
		vector<float> weight1, weight;
		int i, j, edge, number, signal, i1, size, signal1;
		double i0, c0, T1;
		double prob, Psum, Q;

		Q = 1 / (r - 1);
		i0 = 1;
		T1 = 0;
		for (i = 0; i<N; i++)
		{
			T1 = T1 + 1 / pow((i + i0), Q);
		}
		c0 = 1 / T1;
		row.clear();
		bit.clear();
		bit1.clear();
		for (i = 0; i<N; i++)
		{
			bit.push_back(row);
			bit1.push_back(row);

		}
		weight.clear();
		for (i = 1; i <= N; i++)
		{
			weight.push_back(c0 / pow((i + i0 - 1), Q));
		}
		T1 = 0;
		for (i = 0; i<N; i++)
		{
			T1 = T1 + weight[i];
		}
		for (i = 0; i<N; i++)
		{
			weight[i] = weight[i] / T1;
		}

		number = k*N / 2;

		edge = 0;
		while (edge<number)
		{
			nextposition.clear();
			for (j = 0; j<2; j++)
			{
				prob = rand() / float(RAND_MAX);
				// cout<<"prob="<<prob<<endl;
				Psum = 0.0;
				signal = 0;
				size = weight.size();
				for (i1 = 0; i1<size; i1++)
				{
					Psum = Psum + weight[i1];
					if (prob<Psum)
					{
						nextposition.push_back(i1);
						signal = 1;
						break;
					}

				}
				if (signal == 0)
				{
					nextposition.push_back(i1 - 1);
				}
			}

			if (nextposition[0] != nextposition[1])
			{
				signal1 = 0;
				for (intiter = bit1[nextposition[1]].begin(); intiter != bit1[nextposition[1]].end(); intiter++)
				{
					if (*intiter == nextposition[0])
					{
						signal1 = 1;
						break;
					}
				}
				if (signal1 == 0)
				{
					bit1[nextposition[1]].push_back(nextposition[0]);
					bit1[nextposition[0]].push_back(nextposition[1]);
					bit[nextposition[0]].push_back(nextposition[1]);
					edge++;
				}
			}
		}

		vvt::iterator intiter11;
		for (i = 1, intiter11 = bit.begin(); intiter11 != bit.end(); i++, intiter11++)
		{
			row.clear();
			row = *intiter11;
			for (intiter = row.begin(); intiter != row.end(); intiter++)
			{
				/*
				printf("%u %u\n",i,((*intiter)+1));
				ofstream ftest("graph.txt", ios::app);
				ftest<<i<<" "<<((*intiter)+1)<<endl;
				ftest.close();
				*/
				Index* tmp = new Index;
				tmp->x = i;
				tmp->y = ((*intiter) + 1);
				mytest.push_back(*tmp);
				delete tmp;
			}

		}
	//	puts("-1 -1");

	}
	return mytest;
}


class packet {
public:
	int destination; //Ŀ�ĵ�
	int current;   //��ǰλ��
	int downcurrent; //�²�ĵ�ǰλ��
	int downdestination; //�²��Ŀ�ĵ�
	int last; //���ݰ���һ��λ��
	int T; //���ݰ��Ķ���ʱ��
};

int networksize()
{
	int max = 0;
	FILE *p;
	if ((fopen_s(&p,"graph.txt", "rb")) == NULL)
	{
	//	cout << "�Ҳ����ļ�";
	}
	int x, y;
	while (1)
	{
		if (fscanf_s(p, "%d%d", &x, &y,sizeof(int)) == EOF)
			break;
		if (x>max)
		{
			max = x;
		}
		else
		{
		}
		if (y>max)
		{
			max = y;
		}
		else
		{
		}
	}
	return max + 1; //�Ƚϳ������������ڵ�ı�ţ���Ϊ�����Ŵ�0��ʼ���������СΪ�����+1
}
void pathsreen(int from,int to,vector< vector<int> > &path)
{  
	cout << "��" << from << "��" << to << "��·�ɱ�Ϊ" << endl;
	int k;
	k = from;
	if (path[k][to] == -1)
	{
		cout << "����ͨ" << endl;
	}
	else
	{
		if (path[k][to] == to)
		{
			cout << "����ֱ���ʹ�" << endl;
		}
		else
		{

			while (k != to)
			{
				if (path[k][to] == to)
				{

				}
				else
				{
					cout << path[k][to] << endl;
				}
				k = path[k][to];
			}
		}
	}
}
vector< vector<int> >  floyd(vector< vector<int> > &a, int N,double alpha)    //floyd�㷨���·�ɱ�k^alpha
{    
	vector< vector<int> > path(N);  //·�ɱ�
	vector< vector<double> > dis(N);  //�����
	for (int pathstart = 0; pathstart<N; pathstart++)   //��ʼ��·�ɱ���ȫ��-1
	{
		for (int pathdeep = 0; pathdeep<N; pathdeep++)
		{
			path[pathstart].push_back(-1);
		}
	}
	for (int disstart = 0; disstart<N; disstart++)   //��ʼ���������ȫ��MAX
	{
		for (int disdeep = 0; disdeep<N; disdeep++)
		{
			dis[disstart].push_back(M_A_X);
		}
	}
	int have_flag = 0;
	for (int dislong = 0; dislong<N; dislong++)   //��ά�����·�ɱ�
	{
		for (int disdeep = 0; disdeep<N; disdeep++)
		{
			for (int deep = 0; deep < a[dislong].size(); deep++)
			{
				if (a[dislong][deep] == disdeep)
				{
					have_flag = 1;
				}
			}
			if (have_flag == 1)
			{
				dis[dislong][disdeep] = 0;
				path[dislong][disdeep] = disdeep;
			}
			have_flag = 0;
		}
	}
	for (int k = 0; k < N; k++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j<N; j++)
			{
				if (!(dis[i][k] == M_A_X || dis[k][j] == M_A_X) && dis[i][j] >dis[i][k] + dis[k][j]+ pow(a[k].size(),alpha)&&i!=j) 
				{
				    dis[i][j] = dis[i][k] + dis[k][j]+pow(a[k].size(),alpha);
				//	cout << i << " " << j << "���¾���Ϊ" << dis[i][j] <<"����"<<k<< endl;
					path[i][j] = path[i][k];
				}
			}
		}
	}
	return path;
}

vector< vector<int> >  shortest(vector< vector<int> > &a, int N,vector<int> &b)    //floyd�㷨������·��
{
	vector< vector<int> > path(N);  //·�ɱ�
	vector< vector<double> > dis(N);  //�����
	for (int pathstart = 0; pathstart<N; pathstart++)   //��ʼ��·�ɱ���ȫ��-1
	{
		for (int pathdeep = 0; pathdeep<N; pathdeep++)
		{
			path[pathstart].push_back(-1);
		}
	}
	for (int disstart = 0; disstart<N; disstart++)   //��ʼ���������ȫ��MAX
	{
		for (int disdeep = 0; disdeep<N; disdeep++)
		{
			dis[disstart].push_back(M_A_X);
		}
	}
	int have_flag = 0;
	for (int dislong = 0; dislong<N; dislong++)   //��ά�����·�ɱ�
	{
		for (int disdeep = 0; disdeep<N; disdeep++)
		{
			for (int deep = 0; deep < a[dislong].size(); deep++)
			{
				if (a[dislong][deep] == disdeep&&b[dislong]!=1&&b[disdeep]!=1)
				{
					have_flag = 1;
				}
			}
			if (have_flag == 1)
			{
				dis[dislong][disdeep] = 0;
				path[dislong][disdeep] = disdeep;
			}
			have_flag = 0;
		}
	}
	for (int k = 0; k < N; k++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j<N; j++)
			{
				if (!(dis[i][k] == M_A_X || dis[k][j] == M_A_X) && dis[i][j] > dis[i][k] + dis[k][j] + 1) 
				{
					dis[i][j] = dis[i][k] + dis[k][j] + 1;
					//	cout << i << " " << j << "���¾���Ϊ" << dis[i][j] <<"����"<<k<< endl;
					path[i][j] = path[i][k];
				}
			}
		}
	}
	return path;
}




bool compare(int a, int b)
{
	return a>b; //�������У������Ϊreturn a>b����Ϊ����
}





int _tmain(int argc, _TCHAR* argv[])
{
	///////////////////////////////////////////////////////�����ʱ
	clock_t start, finish;
	double totaltime;
	start = clock();
	////////////////////////////////////////////////////////////////////////
	srand((unsigned)time(NULL));   //��ʼ�����������
	int T = 0;//ͳ��������ִ���˶��ٲ�
	double alpha = -3.1; //��ƫ������ߵĲ���
	int graphtime = 0; //ͳ�Ʋ����˶�����ͼ�����ٴζ�����ʵ��
	int  C = 10;// �ڵ�ķ�������
	int N = 1000;  //�����С
	int updiscard = 0;//���ٸ��ϲ����Ϊ�޷��ʹ������
	int downdiscard = 0; //���ٸ��²����Ϊ�޷��ʹﶪ��
	int arrive = 0;   //���ٽڵ��ʹ�
	int clear = 0;
	int bad = 0;
	double lambda = 0.01;   //���ݰ�������
	int overflag = 0; //�Ƿ����ر� 
	int current = 0; //�����ݴ洦��ڵ�ĵ�ǰλ��
	int ranocc = 0;  //0-ռ�ݶ������ֵ ֮��������
	int refresh = 10; //�²�·�ɱ�ĸ���Ƶ��
	int end_flag = 0; //����û�еĽ�����־��
	double gamma = 3.0; //��̬ģ�����ɲ���
	double rand_num = 0; //������ߵ������
	double bia = 0;//��ƫ������ߵĲ���
	int can_reach = 0;// һ���ھ��ܷ񵽴�ı�־
	double energy_begin = 0; //������������ 
	double energy_end = 0;   //����������յ�
	int nextloop = 0; //��һ����ַ
	long double allneighbor_sum = 0;  //���нڵ�֮��
	int same_gamma = 5; //ͬһ��gammaѭ�����ٴ�
	int T_sum;
	vector <int> energy_order;
	//////////////////////////////////////////////////////�������
	int E =6000; //�ڵ��������ֵ 
	long double energysum = 0;   //�ھӽڵ�����֮��
	long double degreesum = 0;  //�ڵ�ĶȺ�
	vector<int> energy;         //��ʾ�²�����ڵ������Ķ�ά����
	vector< vector<int> > up(N); //��ʾ�ϲ�����ṹ�Ķ�ά����
	vector< vector<int> > down(N); //��ʾ�²�����ṹ�Ķ�ά����
	int alpha_id = 0;       //alpha �ı��
	int initial_flag = 0; //all_T�����Ƿ��ʼ��
	vector<int> initial; //�����飬��ʼ����
	vector< vector<int> >  all_T;  //��¼���е�alpha������ʱ��
	vector< vector<int> >  all_arrive;
	vector<double> all_alpha;// ��һ������İ�������
	vector<int> mix;//��ʾ���²�����һһӳ���ϵ������  mix[0]=2 ���ϲ��0�Žڵ��Ӧ�²��2�Žڵ�
	vector<int> remix; //mix���� ��remix[2]=0 �²�2�Ŷ�Ӧ�ϲ��0�Žڵ�
	vector<int> occupy; //��ʾδ��ռ�ݵ��²�ڵ���
	vector <int> ready_attack; //��ʾ׼�����ȥ�����Ľڵ�
	vector <int> attack;
	vector <double> attack_value; //�������²�������Ҫ�Ե�����
	double attack_percentage = 0.1; //��ʾ׼�����й����İٷֱ� 
	int ranatt = 0;
	list<packet> queue;  //�ϲ����ݰ�����
	vector<list<packet>> downqueue;  //�²�ڵ����
	vector<list<packet>> nodequeue;  //�ϲ�ڵ����
	for (int questart = 0; questart < N; questart++)
	{
		nodequeue.push_back(queue);
		downqueue.push_back(queue);
	}
	for (int mixstart = 0; mixstart < N; mixstart++)
	{
		mix.push_back(-1);           //mix��ֵȫ����ʼΪ-1
		remix.push_back(-1);
		energy.push_back(E);        //���нڵ�������ʼΪE 
		occupy.push_back(mixstart);  //�����²�ڵ㶼δ��ռ��
	}
	vector<packet> packet_under_go; //����λ�ڶ���ͷ�ϵ�packet��������ʼ��������飩
	///////////////////////////////////////////////////////////////////�����ʼ��
	while (1)                            ///�ظ�����ʵ��
	{
		if (graphtime == 100)
		{
			break;
		}
		else
		{
		}
		graphtime++;
		gamma = 2.5;
		attack_percentage = 0.0000;
		while (1)
		{
			initial_flag = 0;
			attack_percentage = attack_percentage + 0.00000000000000000000000002;
			if (attack_percentage >= 0.300000001)
			{
				break;
			}
			else
			{
			}
			
			int re = 0;  //ͬһ��gammaѭ�����ٱ�
			while (1)  //ͬһ��gammaҲҪѭ��
			{
				re++;
				if (re > same_gamma)
				{
					break;
				}
			/////////////////////////////////////////////////////////////////////////////////ģ������ ���²��������
			occupy.clear();
			ready_attack.clear();
			attack.clear();
			for (int graphclear = 0; graphclear < N; graphclear++)
			{
				up[graphclear].clear();
				down[graphclear].clear();
				mix[graphclear] = -1;
				remix[graphclear] = -1;
				energy[graphclear] = E;
				occupy.push_back(graphclear);
				ready_attack.push_back(graphclear);
				attack.push_back(0);
			}
			////////////////////////////////////////////////////////////////////////////////// ���²�����ģ��
			vector<Index> uptst = test(N, 6, 2.5);
			for (int i = 0; i < uptst.size(); ++i)
			{
				up[uptst[i].x - 1].push_back(uptst[i].y - 1);
				up[uptst[i].y - 1].push_back(uptst[i].x - 1);
				//	cout << uptst[i].x -1<< "  " << uptst[i].y-1 << endl;
			}
			vector<Index> downtst = test(N, 6, gamma);
			for (int i = 0; i < downtst.size(); ++i)
			{
				down[downtst[i].x - 1].push_back(downtst[i].y - 1);
				down[downtst[i].y - 1].push_back(downtst[i].x - 1);
				//	cout << downtst[i].x -1<< "  " << downtst[i].y-1 << endl;
			}

			///////////////////////////////////////////////////////////////////////////////�����		
			/*
			for (int n = 0; n < N; n++)
			{
				if (mix[n] == -1)
				{
					ranocc = rand() % occupy.size();
					mix[n] = occupy[ranocc];
					remix[mix[n]] = n;
					std::vector<int>::iterator it = occupy.begin() + ranocc;
					occupy.erase(it);
					//cout << "�ϲ�ڵ�" << n << "��Ӧ�²�ڵ�" << mix[n] << endl;
				}
			}
			*/
			/////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////ͬ��
			/*
			int upsizemax = 0;
			int downsizemax = 0;
			int upbig_num = 0;
			int downbig_num = 0;
			vector<int> up_done; //�ϲ��Ѿ���Խ����ĵ�  1Ϊ�� 0Ϊ����
			vector<int> down_done;  //�²��Ѿ���Խ����ĵ�  1Ϊ�� 0Ϊ����
			up_done.clear();
			down_done.clear();
			for (int i = 0; i < N; i++)
			{
			up_done.push_back(0);
			down_done.push_back(0);   //��ʼ��Ϊ0
			}
			for (int n = 0; n < N; n++)
			{
			upsizemax = 0;
			downsizemax = 0;
			upbig_num = -1;
			downbig_num = -1;
			for (int degree_com = 0; degree_com < N; degree_com++)
			{
			if ((up[degree_com].size() >= upsizemax)&&(up_done[degree_com] == 0))
			{
			upsizemax = up[degree_com].size();
			upbig_num = degree_com;
			}
			if ((down[degree_com].size() >= downsizemax)&&(down_done[degree_com] == 0))
			{
			downsizemax = down[degree_com].size();
			downbig_num = degree_com;
			}
			}
			mix[upbig_num] = downbig_num;
			remix[mix[upbig_num]] = upbig_num;
			up_done[upbig_num] = 1;
			down_done[downbig_num] = 1;
			}
			*/
			/*
			for (int n = 0; n < N; n++)
			{
			cout << "�ϲ�ڵ�" << n << "��Ӧ�²�ڵ�" << mix[n] << endl;
			cout << "�²�ڵ�" << n << "��Ӧ�ϲ�ڵ�" << remix[n] << endl;
			system("pause");
			}
			*/
			/////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////����
			
			int upsizemax = 0;
			int downsizemax = 0;
			int upbig_num = 0;
			int downbig_num = 0;
			vector<int> up_done; //�ϲ��Ѿ���Խ����ĵ�  1Ϊ�� 0Ϊ����
			vector<int> down_done;  //�²��Ѿ���Խ����ĵ�  1Ϊ�� 0Ϊ����
			up_done.clear();
			down_done.clear();
			for (int i = 0; i < N; i++)
			{
			up_done.push_back(0);
			down_done.push_back(0);   //��ʼ��Ϊ0
			}
			for (int n = 0; n < N; n++)
			{
			upsizemax = N+1;
			downsizemax = 0;
			upbig_num = -1;
			downbig_num = -1;
			for (int degree_com = 0; degree_com < N; degree_com++)
			{
			if ((up[degree_com].size() <= upsizemax) && (up_done[degree_com] == 0))
			{
			upsizemax = up[degree_com].size();
			upbig_num = degree_com;
			}
			if ((down[degree_com].size() >= downsizemax) && (down_done[degree_com] == 0))
			{
			downsizemax = down[degree_com].size();
			downbig_num = degree_com;
			}
			}
			mix[upbig_num] = downbig_num;
			remix[mix[upbig_num]] = upbig_num;
			up_done[upbig_num] = 1;
			down_done[downbig_num] = 1;
			}
			
			/////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////�����ϲ�·�ɱ��ϲ�·��ֻ����һ��
			for (int rebot = 0; rebot < N; rebot++)      /////////�������
			{
				nodequeue[rebot].clear();
				downqueue[rebot].clear();
			}
			packet_under_go.clear();                     ///����������
			T = 0;//���������
			updiscard = 0;//�����ϲ����Ϊ�޷��ʹ������
			downdiscard = 0; //�²����ݰ�
			current = 0; //�����ݴ洦��ڵ�ĵ�ǰλ��
			alpha = 0.00;
			alpha_id = 0;
			while (1)   //alpha ȡֵ�ı�
			{
				if (initial_flag == 0)
				{
					all_T.push_back(initial);
					all_T[alpha_id].clear();
					all_arrive.push_back(initial);
					all_arrive[alpha_id].clear();
				}
				if (alpha > 1.000005)
				{
					initial_flag = 1;
					break;
				}
				attack_value.clear();
				attack.clear();
				for (int graphclear = 0; graphclear < N; graphclear++)
				{
					attack.push_back(0);
				}
				for (int cal_at_val = 0; cal_at_val < N; cal_at_val++)
				{
					attack_value.push_back(pow(up[cal_at_val].size(), alpha)*pow(down[mix[cal_at_val]].size(), 1 - alpha));
				}
				double max_value = 0;
				int max_value_id = 0;
				for (int att = 0; att < (int)N*attack_percentage; att++)
				{
					max_value = 0;
				    max_value_id = 0;
					for (int i = 0; i < N; i++)
					{
						if (attack_value[i]>max_value)
						{
							max_value = attack_value[i];
							max_value_id = i;
						}					
					}
					attack[max_value_id] = 1;
					//cout << alpha << "�����Ľڵ���Ϊ" << max_value_id << "�����Ϊ" << attack_value[max_value_id] << "�ϲ�Ķ�Ϊ" << up[max_value_id].size() << "�²��Ϊ" << down[mix[max_value_id]].size() << endl;
					attack_value[max_value_id] = 0;
				}
				vector< vector<int> >  uppath = shortest(up, N, attack);      ////////////////////////ѡ��·�ɱ��캯��
				int ret = 0; //ͬһ·�ɱ��ظ����ٴζ���ʵ��
				//////////////////////////////////////////////////////////////////��ձ���
				while (1)
				{
					double true_alpha = 0.00;
					ret++;
					if (ret > 1)
					{
						break;
					}
					for (int rebot = 0; rebot < N; rebot++)      /////////�������,�ڵ���������ֵ
					{
						nodequeue[rebot].clear();
						downqueue[rebot].clear();
						energy[rebot] = E;
					}
					packet_under_go.clear();                     ///����������
					T = 0;//���������
					updiscard = 0;//�����ϲ����Ϊ�޷��ʹ������
					downdiscard = 0; //�²����ݰ�
					current = 0; //�����ݴ洦��ڵ�ĵ�ǰλ��
					end_flag = 0;
					////////////////////////////////////////////////////////////////////////������ʼ
					while (1)
					{
						T++; //���������
						////////////////////////////////////////////////////////////////////////////////////��������
						current = 0; //�����ݴ洦��ڵ�ĵ�ǰλ��
						///////////////////////////////////////////////////////////////////////////////////////				
						///////////////////////////////////////////////////////////////////�������ݰ�
						for (int start_packet = 0; start_packet < (int)(lambda*N*(1-attack_percentage)); start_packet++)  //����lambda*N�����ݰ�
						{
							packet realpacket;
							realpacket.T = 1;
							realpacket.current = rand() % N;  //������ɵ�ǰλ��
							while (attack[realpacket.current] == 1)  //�����ǰλ�ñ�������
							{
								realpacket.current = rand() % N;
							}
							realpacket.destination = realpacket.current;
							realpacket.last = -1;
							realpacket.downcurrent = -1;
							realpacket.downdestination = -1;
							while (realpacket.current == realpacket.destination || attack[realpacket.destination] == 1)  //���Ŀ�ĵر������ˣ��������ɵ�Ŀ�ĵ��뵱ǰλ����ͬ
							{
								realpacket.destination = rand() % N;
							}
							while (uppath[realpacket.current][realpacket.destination]==-1)
							{
								realpacket.current = rand() % N;  //������ɵ�ǰλ��
								while (attack[realpacket.current] == 1)  //�����ǰλ�ñ�������
								{
									realpacket.current = rand() % N;
								}
								realpacket.destination = realpacket.current;
								realpacket.last = -1;
								realpacket.downcurrent = -1;
								realpacket.downdestination = -1;
								while (realpacket.current == realpacket.destination || attack[realpacket.destination] == 1)  //���Ŀ�ĵر������ˣ��������ɵ�Ŀ�ĵ��뵱ǰλ����ͬ
								{
									realpacket.destination = rand() % N;
								}
							}
							nodequeue[realpacket.current].push_back(realpacket);
							//cout<<"���д�С"<<nodequeue[realpacket.current].size() << endl;
						}
						////////////////////////////////////////////////////////////////////////�����ϲ����ݰ�

						for (int search = 0; search < N; search++)
						{
							while (nodequeue[search].size() != 0)//Ҫ��������������ݰ�
							{
								if (uppath[nodequeue[search].front().current][nodequeue[search].front().destination] == -1)
								{
									//cout << "�ϲ�û���ʹ��·�������ݰ�����" << endl;
									updiscard++;
								}
								else
								{
									nodequeue[search].front().downcurrent = mix[search];          //�͵��²�֮ǰ���²�Ķ�Ӧ��ǰλ�øĳ�ӳ��ֵ
									nodequeue[search].front().downdestination = mix[uppath[nodequeue[search].front().current][nodequeue[search].front().destination]];      //�͵��²�֮ǰ���²��Ӧ��Ŀ��λ�øĳ��ϲ�·�ɱ������һ��
									downqueue[mix[search]].push_back(nodequeue[search].front());  //��λ�ڶ���ͷ�ϵ����ݰ��͵��²�
								}
								nodequeue[search].pop_front();
							}
						}
						/*
						for (int search = 0; search < N; search++)
						{
						for (int j = 0; j < C; j++)
						{
						if (nodequeue[search].size() != 0)//Ҫ��������������ݰ�
						{
						packet_under_go.push_back(nodequeue[search].front());  //���ƽ�������У�ԭ�����е�ɾ��
						nodequeue[search].pop_front();
						}
						}
						}
						*/
						packet_under_go.clear();
						////////////////////////////////////////////////////////////���²�����ݰ����д���
						/*
						for (int search = 0; search < N; search++)
						{
						for (int j = 0; j < C; j++)
						{
						if (downqueue[search].size() != 0)//Ҫ��������������ݰ�
						{
						packet_under_go.push_back(downqueue[search].front());  //���ƽ�������У�ԭ�����е�ɾ��
						downqueue[search].pop_front();
						}
						}
						}
						*/
						for (int search = 0; search < N; search++)
						{
							while (downqueue[search].size() != 0)//Ҫ��������������ݰ�
							{
								packet_under_go.push_back(downqueue[search].front());  //���ƽ�������У�ԭ�����е�ɾ��
								downqueue[search].pop_front();
							}
						}
						//////////////////////////////////////////////////////////////////////////////////////////////////////
						for (int process_packet = 0; process_packet < packet_under_go.size(); process_packet++)  //��ʼ��packet undergo���д���
						{
							can_reach = 0;
							if (down[packet_under_go[process_packet].downcurrent].size() == 0||(down[packet_under_go[process_packet].downcurrent].size() == 1 && (attack[down[packet_under_go[process_packet].downcurrent][0]] == 1)))
							{
								//cout << "û���ھӣ����ݰ�����" << endl;
								downdiscard++;
							}
							else
							{
								if (down[packet_under_go[process_packet].downcurrent].size() == 1 && (attack[down[packet_under_go[process_packet].downcurrent][0]]==0))
								{
									energy[packet_under_go[process_packet].downcurrent] = energy[packet_under_go[process_packet].downcurrent] - 1;  //ת���ڵ���ʧ����
									packet_under_go[process_packet].T++;
									packet_under_go[process_packet].last = packet_under_go[process_packet].downcurrent;
									packet_under_go[process_packet].downcurrent = down[packet_under_go[process_packet].downcurrent][0];     //�²㵱ǰλ�ø�ΪΨһ���ھ�
									if (packet_under_go[process_packet].downcurrent == packet_under_go[process_packet].downdestination)   //����²�����ʹ�
									{								
										if (remix[packet_under_go[process_packet].downdestination] == packet_under_go[process_packet].destination) //����²��ʹ�����ϲ�Ŀ�ĵ���ͬ
										{
											arrive++;    //ֱ��ɾ�����ݰ�
											T_sum = T_sum + packet_under_go[process_packet].T;
										}
										else
										{
											packet_under_go[process_packet].current = remix[packet_under_go[process_packet].downdestination];
											nodequeue[remix[packet_under_go[process_packet].downdestination]].push_back(packet_under_go[process_packet]);  //ѹ���ϲ�Ķ���
										}
									}
									else
									{
										downqueue[packet_under_go[process_packet].downcurrent].push_back(packet_under_go[process_packet]);  //ѹ���²�Ķ���
									}
								}
								else
								{
									current = packet_under_go[process_packet].downcurrent; //�ݴ浱ǰλ��
									for (int search_neighbor = 0; search_neighbor < down[current].size(); search_neighbor++) //��ѯ�ھ��Ƿ���ֱ���ʹ�
									{
										if (down[current][search_neighbor] == packet_under_go[process_packet].downdestination)  //����ֱ���ʹ�
										{
											packet_under_go[process_packet].downcurrent = packet_under_go[process_packet].downdestination; //��ΪĿ�ĵ�ַ
											if (remix[packet_under_go[process_packet].downdestination] == packet_under_go[process_packet].destination) //����²��ʹ�����ϲ�Ŀ�ĵ���ͬ
											{
												arrive++;    //ֱ��ɾ�����ݰ�
												T_sum = T_sum + packet_under_go[process_packet].T;
											}
											else
											{
												packet_under_go[process_packet].current = remix[packet_under_go[process_packet].downdestination];
												nodequeue[remix[packet_under_go[process_packet].downdestination]].push_back(packet_under_go[process_packet]);  //ѹ���ϲ�Ķ���
											}
											can_reach = 1;
										}
									}
									if (can_reach == 0)
									{
										for (int allneighbor_count = 0; allneighbor_count < down[current].size(); allneighbor_count++) //ͳ�ƽڵ��ھӵ�������
										{
											if (down[current][allneighbor_count] != packet_under_go[process_packet].last&&attack[down[current][allneighbor_count]]==0)
											{
												allneighbor_sum = allneighbor_sum + pow(down[down[current][allneighbor_count]].size(), true_alpha);
											}
										}
										if (allneighbor_sum==0)
										{
											//cout << "�²������ھӶ���������" << endl;
											if (packet_under_go[process_packet].last==-1)
											{
												downdiscard++;
											}
											else
											{
												packet_under_go[process_packet].T++;
												packet_under_go[process_packet].downcurrent = packet_under_go[process_packet].last;//�޸ĵ�ǰ��ַ
												packet_under_go[process_packet].last = current;
												downqueue[packet_under_go[process_packet].downcurrent].push_back(packet_under_go[process_packet]);  //ѹ���²�Ķ���
												current = -1;
												allneighbor_sum = 0;
												nextloop = 0;
											}
										}
										else
										{
											rand_num = (double)rand() / RAND_MAX;
											for (int deg = 0; deg<down[current].size(); deg++) //ͳ�ƽڵ��ھӵ�����
											{
												if (down[current][deg] != packet_under_go[process_packet].last&&attack[down[current][deg]] == 0)
												{
													energy_begin = energy_end;
													energy_end = energy_end + pow(down[down[current][deg]].size(), true_alpha) / allneighbor_sum;
													//	cout << current << "�ھ���" << down[current][deg] << "��������Ϊ" << energy[down[current][deg]] << "Ϊ����������" << energy_begin << "-" << energy_end << endl;
													if (rand_num>energy_begin&&rand_num < energy_end)
													{
														nextloop = deg; break;
													}
													else
													{
													}
												}
											}
											energy_begin = 0; //��������
											energy_end = 0;
											//	cout << "��һ����ַ" << down[current][nextloop] << "��һ��λ��Ϊ" << packet_under_go[process_packet].last<< endl;
											energy[current] = energy[current] - 1;    //ת����ʧ����
											packet_under_go[process_packet].T++;
											packet_under_go[process_packet].last = current;
											packet_under_go[process_packet].downcurrent = down[current][nextloop];//�޸ĵ�ǰ��ַ				
											downqueue[packet_under_go[process_packet].downcurrent].push_back(packet_under_go[process_packet]);  //ѹ���²�Ķ���
											current = -1;
											energysum = 0;
											degreesum = 0;
											allneighbor_sum = 0;
											nextloop = 0;
										}
									}
								}
							}
						} //packet_under_go ����
						//////////////////////////////////////////////////////////////////////////
						packet_under_go.clear(); //�����������
						//////////////////////////////////////////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////
						/*
						cout << "������е�" << T << "��" << endl;
						cout << "������" << updiscard << "�������ϲ㶪��" << endl;
						cout << "������" << downdiscard << "�������²㶪��" << endl;
						cout <<"������" << arrive << "�����ݰ����ʹ�" << endl;
						*/

						for (int esearch = 0; esearch < N; esearch++)
						{
							if (energy[esearch] <= 0)
							{
								end_flag = 1;
							}
						}
						if (end_flag == 1)
						{
							int energy_min = 10000;
							energy_order.clear();
							double variance_percentage=0.006;
							int energy_sum = 0;
							double variance = 0;
							for (int j = 0; j < N; j++)
							{
								energy_min = 10000;
								for (int i = 0; i < N; i++)
								{
									if (energy[i] < energy_min)
									{
										energy_min = energy[i];
									}
								}
								for (int i = 0; i < N; i++)
								{
									if (energy[i] <= energy_min)
									{
									   energy[i]=6000;
									}
								}
								energy_order.push_back(energy_min);
							}


							for (int i = 0; i < (int)N*variance_percentage; i++)
							{
								energy_sum = energy_order[i] + energy_sum;
							}

							for (int i = 0; i < (int)N*variance_percentage; i++)
							{
								variance = variance+pow(energy_order[i] - energy_sum / ((double)N*variance_percentage), 2);
							}
							cout << "��͵ļ�������" << energy_order[0] << " " << energy_order[1] << " " << energy_order[2] << endl;
							cout << "ƽ������" << energy_sum / ((double)N*variance_percentage) << endl;
							ofstream gtest("finresult.txt", ios::app);
							gtest << alpha << " " << T << " " << arrive << " " << sqrt(variance/((double)N*variance_percentage)) << endl;   ///�ı����
							gtest.close();
							all_T[alpha_id].push_back(T);
							all_arrive[alpha_id]. push_back(arrive);
							break;
						}
						else
						{
						}
					}////������ѭ��
				}///////ͬһ·�ɱ�ѭ��100��
				arrive = 0;
				T_sum = 0;
				updiscard = 0;
				downdiscard = 0;
				if (initial_flag == 0)
				{
					all_alpha.push_back(alpha);
				}
				alpha_id++;
				alpha = alpha + 0.02;
			  }//////////alphaȡֵ�ı�
		    }//////////////һ��gamma��Ҫ��N��
		   int sum_gamma_T = 0; //��ʱ�ܺ�
		   int sum_gamma_T_max = 0; //���ֵ
		   double sum_gamma_rate = 0; //��ʱ�����ܺ�
		   double sum_gamma_rate_max = 0;//�������ֵ
		   int sum_gamma_arrive = 0; //��ʱ�ܺ�
		   int sum_gamma_arrive_max = 0; //���ֵ
		   int max_alpha_id = 0;
		   int max_arrive_id = 0;
		   int max_rate_id = 0;

		   double ave_T = 0;  //ƽ��ֵ
		   double ave_arrive = 0;
		   double ave_rate = 0;
		   int div_count = 0;
		   for (int i = 0; i < all_alpha.size(); i++)
		   {
			   sum_gamma_T = 0;
			   sum_gamma_arrive = 0;
			   sum_gamma_rate = 0;
			   div_count=0;
			   for (int j = 0; j < all_T[i].size(); j++)
			   {
				   sum_gamma_T = sum_gamma_T + all_T[i][j];
				   sum_gamma_arrive = sum_gamma_arrive + all_arrive[i][j];
				   sum_gamma_rate = sum_gamma_rate + all_arrive[i][j]/(all_T[i][j]*lambda*N*(1 - attack_percentage));
				   div_count++;
			   }
			   if (sum_gamma_T>sum_gamma_T_max)
			   {
				   sum_gamma_T_max = sum_gamma_T;
				   max_alpha_id = i;
			   }
			   if (sum_gamma_arrive>sum_gamma_arrive_max)
			   {
				   sum_gamma_arrive_max = sum_gamma_arrive;
				   max_arrive_id = i;
			   }
			   if (sum_gamma_rate>sum_gamma_rate_max)
			   {
				   sum_gamma_rate_max = sum_gamma_rate;
				   max_rate_id = i;
			   }
			   ave_T = (double)sum_gamma_T_max / same_gamma;
			   ave_arrive = (double)sum_gamma_arrive_max / same_gamma;
			   ave_rate = (double)sum_gamma_rate_max / same_gamma;
		   }
		  //cout << "�����ٷֱ�" << attack_percentage << "���²���Ҫ�԰�����" << all_alpha[max_alpha_id] << "��Ӧ�����ʱ��Ϊ" << ave_T << endl;
		   ofstream ftest("gmmaresult.txt", ios::app);
		   ftest << attack_percentage << " " << all_alpha[max_alpha_id] << " " << ave_T << " " << all_alpha[max_arrive_id] << " " << ave_arrive <<" "<<all_alpha[max_rate_id]<<" "<<ave_rate<< endl;   ///�ı����
		   ftest.close();
		}//gamma ȡֵ�ı�
		////////////////////////////////////////////////////////////////////////��ʱ����
		finish = clock();
		totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "\n�˳�������һ�ζ���ʵ���ʱ��Ϊ" << totaltime << "��" << endl;
	}
}



		

