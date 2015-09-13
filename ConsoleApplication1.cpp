// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
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
#define M_A_X 9999999999999999999  //MAX定义为30万
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
	int destination; //目的地
	int current;   //当前位置
	int downcurrent; //下层的当前位置
	int downdestination; //下层的目的地
	int last; //数据包上一跳位置
	int T; //数据包的逗留时间
};

int networksize()
{
	int max = 0;
	FILE *p;
	if ((fopen_s(&p,"graph.txt", "rb")) == NULL)
	{
	//	cout << "找不到文件";
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
	return max + 1; //比较出网络里面最大节点的编号，因为网络编号从0开始的则网络大小为最大编号+1
}
void pathsreen(int from,int to,vector< vector<int> > &path)
{  
	cout << "从" << from << "到" << to << "的路由表为" << endl;
	int k;
	k = from;
	if (path[k][to] == -1)
	{
		cout << "不连通" << endl;
	}
	else
	{
		if (path[k][to] == to)
		{
			cout << "可以直接送达" << endl;
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
vector< vector<int> >  floyd(vector< vector<int> > &a, int N,double alpha)    //floyd算法输出路由表，k^alpha
{    
	vector< vector<int> > path(N);  //路由表
	vector< vector<double> > dis(N);  //距离表
	for (int pathstart = 0; pathstart<N; pathstart++)   //初始化路由表组全是-1
	{
		for (int pathdeep = 0; pathdeep<N; pathdeep++)
		{
			path[pathstart].push_back(-1);
		}
	}
	for (int disstart = 0; disstart<N; disstart++)   //初始化距离表组全是MAX
	{
		for (int disdeep = 0; disdeep<N; disdeep++)
		{
			dis[disstart].push_back(M_A_X);
		}
	}
	int have_flag = 0;
	for (int dislong = 0; dislong<N; dislong++)   //二维数组存路由表
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
				//	cout << i << " " << j << "的新距离为" << dis[i][j] <<"经过"<<k<< endl;
					path[i][j] = path[i][k];
				}
			}
		}
	}
	return path;
}

vector< vector<int> >  shortest(vector< vector<int> > &a, int N,vector<int> &b)    //floyd算法输出最短路径
{
	vector< vector<int> > path(N);  //路由表
	vector< vector<double> > dis(N);  //距离表
	for (int pathstart = 0; pathstart<N; pathstart++)   //初始化路由表组全是-1
	{
		for (int pathdeep = 0; pathdeep<N; pathdeep++)
		{
			path[pathstart].push_back(-1);
		}
	}
	for (int disstart = 0; disstart<N; disstart++)   //初始化距离表组全是MAX
	{
		for (int disdeep = 0; disdeep<N; disdeep++)
		{
			dis[disstart].push_back(M_A_X);
		}
	}
	int have_flag = 0;
	for (int dislong = 0; dislong<N; dislong++)   //二维数组存路由表
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
					//	cout << i << " " << j << "的新距离为" << dis[i][j] <<"经过"<<k<< endl;
					path[i][j] = path[i][k];
				}
			}
		}
	}
	return path;
}




bool compare(int a, int b)
{
	return a>b; //升序排列，如果改为return a>b，则为降序
}





int _tmain(int argc, _TCHAR* argv[])
{
	///////////////////////////////////////////////////////程序计时
	clock_t start, finish;
	double totaltime;
	start = clock();
	////////////////////////////////////////////////////////////////////////
	srand((unsigned)time(NULL));   //初始化随机数种子
	int T = 0;//统计主程序执行了多少步
	double alpha = -3.1; //有偏随机行走的参数
	int graphtime = 0; //统计产生了多少张图，多少次独立的实验
	int  C = 10;// 节点的发送能力
	int N = 1000;  //网络大小
	int updiscard = 0;//多少个上层包因为无法送达而丢弃
	int downdiscard = 0; //多少个下层包因为无法送达丢弃
	int arrive = 0;   //多少节点送达
	int clear = 0;
	int bad = 0;
	double lambda = 0.01;   //数据包生成率
	int overflag = 0; //是否有重边 
	int current = 0; //用来暂存处理节点的当前位置
	int ranocc = 0;  //0-占据队列最大值 之间的随机数
	int refresh = 10; //下层路由表的更新频率
	int end_flag = 0; //能量没有的结束标志符
	double gamma = 3.0; //静态模型生成参数
	double rand_num = 0; //随机行走的随机数
	double bia = 0;//有偏随机行走的参数
	int can_reach = 0;// 一跳邻居能否到达的标志
	double energy_begin = 0; //能量分配的起点 
	double energy_end = 0;   //能量分配的终点
	int nextloop = 0; //下一跳地址
	long double allneighbor_sum = 0;  //所有节点之和
	int same_gamma = 5; //同一个gamma循环多少次
	int T_sum;
	vector <int> energy_order;
	//////////////////////////////////////////////////////网络参数
	int E =6000; //节点的能量初值 
	long double energysum = 0;   //邻居节点能量之合
	long double degreesum = 0;  //节点的度和
	vector<int> energy;         //表示下层网络节点能量的二维数组
	vector< vector<int> > up(N); //表示上层网络结构的二维数组
	vector< vector<int> > down(N); //表示下层网络结构的二维数组
	int alpha_id = 0;       //alpha 的编号
	int initial_flag = 0; //all_T数组是否初始化
	vector<int> initial; //空数组，初始化用
	vector< vector<int> >  all_T;  //记录所有的alpha下生存时间
	vector< vector<int> >  all_arrive;
	vector<double> all_alpha;// 上一个数组的伴生数组
	vector<int> mix;//表示上下层网络一一映射关系的数组  mix[0]=2 即上层的0号节点对应下层的2号节点
	vector<int> remix; //mix的逆 即remix[2]=0 下层2号对应上层的0号节点
	vector<int> occupy; //表示未被占据的下层节点编号
	vector <int> ready_attack; //表示准备抽调去攻击的节点
	vector <int> attack;
	vector <double> attack_value; //度量上下层网络重要性的数组
	double attack_percentage = 0.1; //表示准备进行攻击的百分比 
	int ranatt = 0;
	list<packet> queue;  //上层数据包队列
	vector<list<packet>> downqueue;  //下层节点队列
	vector<list<packet>> nodequeue;  //上层节点队列
	for (int questart = 0; questart < N; questart++)
	{
		nodequeue.push_back(queue);
		downqueue.push_back(queue);
	}
	for (int mixstart = 0; mixstart < N; mixstart++)
	{
		mix.push_back(-1);           //mix的值全部初始为-1
		remix.push_back(-1);
		energy.push_back(E);        //所有节点能量初始为E 
		occupy.push_back(mixstart);  //所有下层节点都未被占据
	}
	vector<packet> packet_under_go; //所有位于队列头上的packet（即将开始处理的数组）
	///////////////////////////////////////////////////////////////////网络初始化
	while (1)                            ///重复进行实验
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
			
			int re = 0;  //同一个gamma循环多少遍
			while (1)  //同一个gamma也要循环
			{
				re++;
				if (re > same_gamma)
				{
					break;
				}
			/////////////////////////////////////////////////////////////////////////////////模型生成 上下层网络清空
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
			////////////////////////////////////////////////////////////////////////////////// 上下层网络模型
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

			///////////////////////////////////////////////////////////////////////////////随机配		
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
					//cout << "上层节点" << n << "对应下层节点" << mix[n] << endl;
				}
			}
			*/
			/////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////同配
			/*
			int upsizemax = 0;
			int downsizemax = 0;
			int upbig_num = 0;
			int downbig_num = 0;
			vector<int> up_done; //上层已经配对结束的点  1为是 0为不是
			vector<int> down_done;  //下层已经配对结束的点  1为是 0为不是
			up_done.clear();
			down_done.clear();
			for (int i = 0; i < N; i++)
			{
			up_done.push_back(0);
			down_done.push_back(0);   //初始化为0
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
			cout << "上层节点" << n << "对应下层节点" << mix[n] << endl;
			cout << "下层节点" << n << "对应上层节点" << remix[n] << endl;
			system("pause");
			}
			*/
			/////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////异配
			
			int upsizemax = 0;
			int downsizemax = 0;
			int upbig_num = 0;
			int downbig_num = 0;
			vector<int> up_done; //上层已经配对结束的点  1为是 0为不是
			vector<int> down_done;  //下层已经配对结束的点  1为是 0为不是
			up_done.clear();
			down_done.clear();
			for (int i = 0; i < N; i++)
			{
			up_done.push_back(0);
			down_done.push_back(0);   //初始化为0
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
			///////////////////////////////////////////////////////////////////构造上层路由表，上层路由只构造一次
			for (int rebot = 0; rebot < N; rebot++)      /////////队列清空
			{
				nodequeue[rebot].clear();
				downqueue[rebot].clear();
			}
			packet_under_go.clear();                     ///处理队列清空
			T = 0;//主程序计数
			updiscard = 0;//多少上层包因为无法送达而丢弃
			downdiscard = 0; //下层数据包
			current = 0; //用来暂存处理节点的当前位置
			alpha = 0.00;
			alpha_id = 0;
			while (1)   //alpha 取值改变
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
					//cout << alpha << "攻击的节点编号为" << max_value_id << "其度量为" << attack_value[max_value_id] << "上层的度为" << up[max_value_id].size() << "下层度为" << down[mix[max_value_id]].size() << endl;
					attack_value[max_value_id] = 0;
				}
				vector< vector<int> >  uppath = shortest(up, N, attack);      ////////////////////////选好路由表构造函数
				int ret = 0; //同一路由表重复多少次独立实验
				//////////////////////////////////////////////////////////////////清空变量
				while (1)
				{
					double true_alpha = 0.00;
					ret++;
					if (ret > 1)
					{
						break;
					}
					for (int rebot = 0; rebot < N; rebot++)      /////////队列清空,节点能量赋初值
					{
						nodequeue[rebot].clear();
						downqueue[rebot].clear();
						energy[rebot] = E;
					}
					packet_under_go.clear();                     ///处理队列清空
					T = 0;//主程序计数
					updiscard = 0;//多少上层包因为无法送达而丢弃
					downdiscard = 0; //下层数据包
					current = 0; //用来暂存处理节点的当前位置
					end_flag = 0;
					////////////////////////////////////////////////////////////////////////主程序开始
					while (1)
					{
						T++; //主程序计数
						////////////////////////////////////////////////////////////////////////////////////变量清零
						current = 0; //用来暂存处理节点的当前位置
						///////////////////////////////////////////////////////////////////////////////////////				
						///////////////////////////////////////////////////////////////////生成数据包
						for (int start_packet = 0; start_packet < (int)(lambda*N*(1-attack_percentage)); start_packet++)  //生成lambda*N个数据包
						{
							packet realpacket;
							realpacket.T = 1;
							realpacket.current = rand() % N;  //随机生成当前位置
							while (attack[realpacket.current] == 1)  //如果当前位置被攻击了
							{
								realpacket.current = rand() % N;
							}
							realpacket.destination = realpacket.current;
							realpacket.last = -1;
							realpacket.downcurrent = -1;
							realpacket.downdestination = -1;
							while (realpacket.current == realpacket.destination || attack[realpacket.destination] == 1)  //如果目的地被攻击了，或是生成的目的地与当前位置相同
							{
								realpacket.destination = rand() % N;
							}
							while (uppath[realpacket.current][realpacket.destination]==-1)
							{
								realpacket.current = rand() % N;  //随机生成当前位置
								while (attack[realpacket.current] == 1)  //如果当前位置被攻击了
								{
									realpacket.current = rand() % N;
								}
								realpacket.destination = realpacket.current;
								realpacket.last = -1;
								realpacket.downcurrent = -1;
								realpacket.downdestination = -1;
								while (realpacket.current == realpacket.destination || attack[realpacket.destination] == 1)  //如果目的地被攻击了，或是生成的目的地与当前位置相同
								{
									realpacket.destination = rand() % N;
								}
							}
							nodequeue[realpacket.current].push_back(realpacket);
							//cout<<"队列大小"<<nodequeue[realpacket.current].size() << endl;
						}
						////////////////////////////////////////////////////////////////////////处理上层数据包

						for (int search = 0; search < N; search++)
						{
							while (nodequeue[search].size() != 0)//要求队列里面有数据包
							{
								if (uppath[nodequeue[search].front().current][nodequeue[search].front().destination] == -1)
								{
									//cout << "上层没有送达的路径，数据包丢弃" << endl;
									updiscard++;
								}
								else
								{
									nodequeue[search].front().downcurrent = mix[search];          //送到下层之前将下层的对应当前位置改成映射值
									nodequeue[search].front().downdestination = mix[uppath[nodequeue[search].front().current][nodequeue[search].front().destination]];      //送到下层之前将下层对应的目的位置改成上层路由表里的下一跳
									downqueue[mix[search]].push_back(nodequeue[search].front());  //将位于队列头上的数据包送到下层
								}
								nodequeue[search].pop_front();
							}
						}
						/*
						for (int search = 0; search < N; search++)
						{
						for (int j = 0; j < C; j++)
						{
						if (nodequeue[search].size() != 0)//要求队列里面有数据包
						{
						packet_under_go.push_back(nodequeue[search].front());  //复制进处理队列，原来队列的删掉
						nodequeue[search].pop_front();
						}
						}
						}
						*/
						packet_under_go.clear();
						////////////////////////////////////////////////////////////对下层的数据包进行处理
						/*
						for (int search = 0; search < N; search++)
						{
						for (int j = 0; j < C; j++)
						{
						if (downqueue[search].size() != 0)//要求队列里面有数据包
						{
						packet_under_go.push_back(downqueue[search].front());  //复制进处理队列，原来队列的删掉
						downqueue[search].pop_front();
						}
						}
						}
						*/
						for (int search = 0; search < N; search++)
						{
							while (downqueue[search].size() != 0)//要求队列里面有数据包
							{
								packet_under_go.push_back(downqueue[search].front());  //复制进处理队列，原来队列的删掉
								downqueue[search].pop_front();
							}
						}
						//////////////////////////////////////////////////////////////////////////////////////////////////////
						for (int process_packet = 0; process_packet < packet_under_go.size(); process_packet++)  //开始对packet undergo进行处理
						{
							can_reach = 0;
							if (down[packet_under_go[process_packet].downcurrent].size() == 0||(down[packet_under_go[process_packet].downcurrent].size() == 1 && (attack[down[packet_under_go[process_packet].downcurrent][0]] == 1)))
							{
								//cout << "没有邻居，数据包丢弃" << endl;
								downdiscard++;
							}
							else
							{
								if (down[packet_under_go[process_packet].downcurrent].size() == 1 && (attack[down[packet_under_go[process_packet].downcurrent][0]]==0))
								{
									energy[packet_under_go[process_packet].downcurrent] = energy[packet_under_go[process_packet].downcurrent] - 1;  //转发节点损失能量
									packet_under_go[process_packet].T++;
									packet_under_go[process_packet].last = packet_under_go[process_packet].downcurrent;
									packet_under_go[process_packet].downcurrent = down[packet_under_go[process_packet].downcurrent][0];     //下层当前位置改为唯一的邻居
									if (packet_under_go[process_packet].downcurrent == packet_under_go[process_packet].downdestination)   //如果下层可以送达
									{								
										if (remix[packet_under_go[process_packet].downdestination] == packet_under_go[process_packet].destination) //如果下层送达后与上层目的地相同
										{
											arrive++;    //直接删除数据包
											T_sum = T_sum + packet_under_go[process_packet].T;
										}
										else
										{
											packet_under_go[process_packet].current = remix[packet_under_go[process_packet].downdestination];
											nodequeue[remix[packet_under_go[process_packet].downdestination]].push_back(packet_under_go[process_packet]);  //压入上层的队列
										}
									}
									else
									{
										downqueue[packet_under_go[process_packet].downcurrent].push_back(packet_under_go[process_packet]);  //压入下层的队列
									}
								}
								else
								{
									current = packet_under_go[process_packet].downcurrent; //暂存当前位置
									for (int search_neighbor = 0; search_neighbor < down[current].size(); search_neighbor++) //查询邻居是否能直接送达
									{
										if (down[current][search_neighbor] == packet_under_go[process_packet].downdestination)  //可以直接送达
										{
											packet_under_go[process_packet].downcurrent = packet_under_go[process_packet].downdestination; //改为目的地址
											if (remix[packet_under_go[process_packet].downdestination] == packet_under_go[process_packet].destination) //如果下层送达后与上层目的地相同
											{
												arrive++;    //直接删除数据包
												T_sum = T_sum + packet_under_go[process_packet].T;
											}
											else
											{
												packet_under_go[process_packet].current = remix[packet_under_go[process_packet].downdestination];
												nodequeue[remix[packet_under_go[process_packet].downdestination]].push_back(packet_under_go[process_packet]);  //压入上层的队列
											}
											can_reach = 1;
										}
									}
									if (can_reach == 0)
									{
										for (int allneighbor_count = 0; allneighbor_count < down[current].size(); allneighbor_count++) //统计节点邻居的能量和
										{
											if (down[current][allneighbor_count] != packet_under_go[process_packet].last&&attack[down[current][allneighbor_count]]==0)
											{
												allneighbor_sum = allneighbor_sum + pow(down[down[current][allneighbor_count]].size(), true_alpha);
											}
										}
										if (allneighbor_sum==0)
										{
											//cout << "下层所有邻居都被攻击了" << endl;
											if (packet_under_go[process_packet].last==-1)
											{
												downdiscard++;
											}
											else
											{
												packet_under_go[process_packet].T++;
												packet_under_go[process_packet].downcurrent = packet_under_go[process_packet].last;//修改当前地址
												packet_under_go[process_packet].last = current;
												downqueue[packet_under_go[process_packet].downcurrent].push_back(packet_under_go[process_packet]);  //压入下层的队列
												current = -1;
												allneighbor_sum = 0;
												nextloop = 0;
											}
										}
										else
										{
											rand_num = (double)rand() / RAND_MAX;
											for (int deg = 0; deg<down[current].size(); deg++) //统计节点邻居的能量
											{
												if (down[current][deg] != packet_under_go[process_packet].last&&attack[down[current][deg]] == 0)
												{
													energy_begin = energy_end;
													energy_end = energy_end + pow(down[down[current][deg]].size(), true_alpha) / allneighbor_sum;
													//	cout << current << "邻居是" << down[current][deg] << "它的能量为" << energy[down[current][deg]] << "为其分配的区间" << energy_begin << "-" << energy_end << endl;
													if (rand_num>energy_begin&&rand_num < energy_end)
													{
														nextloop = deg; break;
													}
													else
													{
													}
												}
											}
											energy_begin = 0; //区间清零
											energy_end = 0;
											//	cout << "下一跳地址" << down[current][nextloop] << "上一次位置为" << packet_under_go[process_packet].last<< endl;
											energy[current] = energy[current] - 1;    //转发损失能量
											packet_under_go[process_packet].T++;
											packet_under_go[process_packet].last = current;
											packet_under_go[process_packet].downcurrent = down[current][nextloop];//修改当前地址				
											downqueue[packet_under_go[process_packet].downcurrent].push_back(packet_under_go[process_packet]);  //压入下层的队列
											current = -1;
											energysum = 0;
											degreesum = 0;
											allneighbor_sum = 0;
											nextloop = 0;
										}
									}
								}
							}
						} //packet_under_go 处理
						//////////////////////////////////////////////////////////////////////////
						packet_under_go.clear(); //处理队列清零
						//////////////////////////////////////////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////
						/*
						cout << "程序进行到" << T << "步" << endl;
						cout << "本次有" << updiscard << "个包被上层丢弃" << endl;
						cout << "本次有" << downdiscard << "个包被下层丢弃" << endl;
						cout <<"本次有" << arrive << "个数据包被送达" << endl;
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
							cout << "最低的几个能量" << energy_order[0] << " " << energy_order[1] << " " << energy_order[2] << endl;
							cout << "平均能量" << energy_sum / ((double)N*variance_percentage) << endl;
							ofstream gtest("finresult.txt", ios::app);
							gtest << alpha << " " << T << " " << arrive << " " << sqrt(variance/((double)N*variance_percentage)) << endl;   ///文本输出
							gtest.close();
							all_T[alpha_id].push_back(T);
							all_arrive[alpha_id]. push_back(arrive);
							break;
						}
						else
						{
						}
					}////主程序循环
				}///////同一路由表循环100次
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
			  }//////////alpha取值改变
		    }//////////////一个gamma下要算N遍
		   int sum_gamma_T = 0; //暂时总和
		   int sum_gamma_T_max = 0; //最大值
		   double sum_gamma_rate = 0; //暂时比率总和
		   double sum_gamma_rate_max = 0;//比率最大值
		   int sum_gamma_arrive = 0; //暂时总和
		   int sum_gamma_arrive_max = 0; //最大值
		   int max_alpha_id = 0;
		   int max_arrive_id = 0;
		   int max_rate_id = 0;

		   double ave_T = 0;  //平均值
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
		  //cout << "攻击百分比" << attack_percentage << "上下层重要性阿尔法" << all_alpha[max_alpha_id] << "对应的最大时间为" << ave_T << endl;
		   ofstream ftest("gmmaresult.txt", ios::app);
		   ftest << attack_percentage << " " << all_alpha[max_alpha_id] << " " << ave_T << " " << all_alpha[max_arrive_id] << " " << ave_arrive <<" "<<all_alpha[max_rate_id]<<" "<<ave_rate<< endl;   ///文本输出
		   ftest.close();
		}//gamma 取值改变
		////////////////////////////////////////////////////////////////////////计时结束
		finish = clock();
		totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "\n此程序运行一次独立实验的时间为" << totaltime << "秒" << endl;
	}
}



		

