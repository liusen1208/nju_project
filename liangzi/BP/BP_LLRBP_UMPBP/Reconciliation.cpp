#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <climits>
#include <cfloat>
#include <vector>
#define n 10000                //H的规模：m*n       //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
#define m 2000                                      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
using namespace std;

/*
变量节点(variance node),校验节点（check node）分别简称为VN,CN ；
校验矩阵，用H表示；   Alice/Bob筛后密钥 ：X/Y；    Alice校验子：Z
*/

//注：“--------------”表示有改进的地方 

//边结构体 
struct Edge {
	int N;        //CN/VN 
	double Lq;    //L(qij) 
	double Lr;    //L(rji)
	int numOfRL;  //Hr[][]（Hl[][]）该元素和 Hl[][numOfRL] (Hr[][numOfRL])相等; Hr中该1是Hl中的对应列的第几个1
};
vector<Edge> Hr[m];         //H的行优先存储
vector<Edge> Hl[n];         //H的列优先存储 
int Dv[n] = { 0 };          //VN度分布 
int Dc[m] = { 0 };          //CN度分布
int Y[n];                   //Y 
int Z[m];                   //Z 
int X[n];                   //X,测试所用，之后删除 


/* 用来确定Hr和Hl中Edge结构体中的 numOfRL值 */    
void calNumOfRL(){                      
	int i,j,k,vn,cn,tmp;
	for(i=0;i<n;i++){      //遍历VN
		cn=Dv[i];                //该VN度 (Dv已经赋初值)
		for(k=0;k<cn;k++){       //遍历与该VN相连的CN 
			j=Hl[i][k].N;              //与VN i相连的CN j 
			tmp = 0;
			while(Hr[j][tmp].N != i) tmp++;  //找出VN i为第几个与CN j相连VN 
			Hl[i][k].numOfRL=tmp;      //求出numOfRL值           
		} 
	}
	cout<<"    变量节点中间量计算完毕！"<<endl; 
	for(j=0;j<m;j++){      //遍历CN
		vn=Dc[j];                //该CN度 (Dc已经赋初值)
		for(k=0;k<vn;k++){       //遍历与该CN相连的VN 
			i=Hr[j][k].N;              //与CN j相连的VN i
			tmp = 0;
			while(Hl[i][tmp].N != j) tmp++;  //找出CN j为第几个与VN i相连CN 
			Hr[j][k].numOfRL=tmp;      //求出numOfRL值           
		} 
	}
	cout<<"    校验节点中间量计算完毕！"<<endl; 
}

/*从文本中读取H，并初始化 Dc 和 Dv*/
void readH(const char * address) {
	ifstream ifile;
	ifile.open(address);
	char buff[50];
	int row, column;
	Edge ed;
	
	cout<<"1.开始从 \""<<address<<"\" 读LDPC码！"<<endl;
	cout<<"    读取中..."<<endl; 
	
	while (ifile.getline(buff, 50, ' ')) {//读出行 
		row = atoi(buff);                 	  //转成整型 
		ifile.getline(buff, 50, '\n');        //读出列 
		column = atoi(buff);                  //转成整型
		ed.N = column; Hr[row].push_back(ed); //放到Hr,Hl 
		ed.N = row; Hl[column].push_back(ed);
		Dc[row]++;                            //该CN的度+1
		Dv[column]++;                         //该VN的度+1
	}
	ifile.close();
	
	cout<<"  LDPC码读取结束！"<<endl<<endl;
}

//测试X和Y有多少不同 
int Test() {  //另一方的最终密钥，长度
	int sum = 0, i;
	for (i = 0;i < n;i++) {
		if (X[i] != Y[i])
			sum++;
	}
	return sum;
}



/*LLR-BP算法，返回协商是否正常结束(正常返回迭代次数，不正常返回-1) <要保证 Dv 和 Dc已经初始化> */
int LLRBP(double e,int iteraM) {   //误码率，Alice校验子，迭代次数 
    //---------------变量中，去除了  tagV[n] 和 tmp00，增加了 tmp2,tmp3,tmp4,tmp5,tmp6------------------------------------ 
	int i, j, k, iteraN = 0, tmpe,tmp2,tmp4;         
	double tmp0, tmp1,tmp3,tmp5,tmp6;            
	bool isEqual;                 //判断H*Y是否等于Z 
	double Lp[n];                 //L(Pi)
	double LQ[n];                 //L(Qi)
	int vn, cn;                   //暂存的VN/CN个数 
	
	//初始化L(Pi)
	double a0 = log((1 - e) / e), a1 = log(e / (1 - e));           //L(Pi)两种取值 
	for (i = 0;i<n;i++)                                            //对每个VN i 
		Lp[i] = (Y[i] == 0) ? a0 : a1;	                           //初始化L(Pi)
		
		
	//初始化L(qij)
	for (j = 0;j<m;j++) {                                          //遍历每个CN
		vn = Dc[j];                                                //与该CN相连VN数(即CN的度)   
		for (i = 0;i<vn;i++) {                                     //遍历与该CN相连VN 
			tmpe = Hr[j][i].N;                                     //VN
			Hr[j][i].Lq = Lp[tmpe];                                //L(qij)	
			Hl[tmpe][Hr[j][i].numOfRL].Lq = Lp[tmpe];              //L(qij)
		}
	}




	
	
	//处理CN/VN消息，译码判决，停止判断 
	do {
		iteraN++;                                                  //迭代次数加1 
		
		//处理CN消息 
		for (j = 0;j<m;j++) {                                      //遍历每个CN
			vn = Dc[j];                                            //与该CN相连VN数  			 
			tmp0 = ((Z[j]==0) ? 1.0 : -1.0);                       //符号函数 

			k=0; 
			for(tmp3=1.0,i=0;i<vn;i++){                            //遍历这些VN，计算它们的双曲正切之积 
				if(Hr[j][i].Lq!=0){                      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					tmp3=tmp3*tanh(0.5*Hr[j][i].Lq);
				}else{
					k++;
				}
			}
			for (i = 0;i<vn;i++) {                                 //遍历这些VN 
				if(k==0){                                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					tmp1=tmp3/tanh(0.5*Hr[j][i].Lq);	                                    //计算公式中两乘积	
				}else{
					if(k==1){
						if(Hr[j][i].Lq==0){
							tmp1=tmp3;
						}else{
							tmp1=0;
						}
					}else{
						tmp1=0;
					}
				}
				Hr[j][i].Lr = tmp0*log((1.0 + tmp1) / (1.0 - tmp1));//L(r j_tmpe)
				tmpe = Hr[j][i].N;                                  //该VN 
				Hl[tmpe][Hr[j][i].numOfRL].Lr = Hr[j][i].Lr;        //L(r j_tmpe)
			}
		}
	
		//处理VN消息，译码判决 	
		for (i = 0;i<n;i++) {           //遍历每个VN
				
			//译码判决该VN 
			cn = Dv[i];                       //与该VN相连CN数   
			for (tmp0 = 0.0,j = 0;j<cn;j++) { //遍历这些CN           
				tmp0 = tmp0+Hl[i][j].Lr;
		    }
		    LQ[i] = Lp[i] + tmp0;         //L(Qi)	
			if (LQ[i]>=0)	Y[i]=0;	
			else            Y[i]=1;
						
			//处理该VN消息
			for(j=0;j<cn;j++){                  //遍历这些CN  
				Hl[i][j].Lq=LQ[i]-Hl[i][j].Lr;    //L(q i_tmpe)
				tmpe = Hl[i][j].N;                      //该CN 
				Hr[tmpe][Hl[i][j].numOfRL].Lq = Hl[i][j].Lq;           //L(q i_tmpe) 
			}   		                             	                             
		}
		
		//停止判断
		isEqual = true;               //设H*Y=Z 
		for (j = 0;j<m;j++) {           //计算H*Y 
			vn=Dc[j];
						
			for (k=0,i = 0;i<vn;i++){       //遍历与 CN j相连的所有 VN 
				tmpe=Hr[j][i].N;                  //对VN tmpe
				if(Y[tmpe]==1) k++;               //H第j行与Y对应位相乘为 1 的个数加1 
			}	
			k=((k%2)==1)?1:0;               //H第j行与Y相乘(为得到0/1结果，*/+分别用&&/^代替)
							
			if (k != Z[j]) {                  //H*Y!=Z 
				isEqual = false;                 //标记H*Y!=Z 
				break;                         //停止剩余计算 
			}
		}
	} while (iteraN<iteraM&&!isEqual);  //未达最大迭代，且H*Y!=Z 
	
	if(!isEqual){
		iteraN++;
	}
	return iteraN;       //返回迭代次数 
}



//获取Alice长度为m校验子Z
void getZ(){
	int i,j,k,vn;
	//------------------------------------------------------
	for (j = 0;j<m;j++) {           //计算H*X		
		vn=Dc[j]; k=0;
		for(i=0;i<vn;i++){          //遍历与 CN j相连的所有 VN 
			if(X[Hr[j][i].N]==1) k++;     //H第j行与Y对应位相乘为 1 的个数加1 
		}	
		Z[j]=((k%2)==1)?1:0;        //H第j行与Y相乘(为得到0/1结果，*/+分别用&&/^代替)	
	}	
} 

//模拟X
void randomX(){
	int i;
	for(i=0;i<n;i++){
		X[i]=rand() % 2;
	}
} 

//由X模拟出Y 
void randomY(double e){
	int errorN=round(n*e);
	int i,k;
	int mark[n];
	for(i=0;i<n;i++){
		Y[i]=X[i];
		mark[i]=0;
	}
	for(i=0;i<errorN;i++){  //调控错误数   
		k=rand()%n;
		if(mark[k]==1){
			i--;
			continue;
		}
		Y[k]=(Y[k]==0)?1:0;
		mark[k]=1;
	}
}

//由码率和错误率，计算协商效率f
double calF(double e){
	double H=-e*(log(e)/log(2))-(1-e)*(log(1-e)/log(2));
	return m*1.0/n/H;
} 

//归纳出一个矩阵的变量节点度分布系数
void concludeDegreeCoefficient(){
	
	int DegreeNum[1000];        //度数i的点数；设度数不超过999
	memset(DegreeNum,0,1000*sizeof(int));  //初始化
	
	for(int i=0;i<n;i++){
		DegreeNum[Dv[i]]++;     //该度数的点数加 1 
	} 
	
	int Sum=0;
	for(int i=0;i<1000;i++){    //计算总边数 
		Sum+=(i*DegreeNum[i]);
	} 
	
	double d=0.0,tmp; 
	for(int i=0;i<1000;i++){    //输出系数，和相对应的度 
		if(DegreeNum[i]!=0){
			tmp=(i*DegreeNum[i]*1.0)/Sum;
			d+=tmp;
			cout<<tmp<<"     "<<i<<endl;
		}
	}
	cout<<d<<endl;
} 



//产生指定条数、指定误码率的密钥X，Y，并写进两个文件中
void createXY(int num,double e){
	int k,kk,i;
	int errorN=round(n*e);
	int mark[n];
	ofstream ofile1,ofile2;
	ofile1.open("X.txt");
	ofile2.open("Y.txt");
	
	for(kk=0;kk<num;kk++){               //产生num条密钥 
		for(i=0;i<(n-1);i++){               //产生密钥的前(n-1)位 
			X[i]=rand() % 2;                      //产生X，写入文件 
			ofile1<<X[i]<<" ";
		}
		X[n-1]=rand() % 2;                  //产生密钥的最后一位            
		ofile1<<X[n-1]<<endl;
		
		for(i=0;i<n;i++){
			Y[i]=X[i];
			mark[i]=0;
		}
		for(i=0;i<errorN;i++){          //调控错误数    
			k=rand()%n;
			if(mark[k]==1){
				i--;
				continue;
			}
			Y[k]=(Y[k]==0)?1:0;
			mark[k]=1;
		}
		for(i=0;i<(n-1);i++){        //Y写入文件 
			ofile2<<Y[i]<<" ";
		} 
		ofile2<<Y[n-1]<<endl;
		
	}
	
	ofile1.close();
	ofile2.close();
} 


//从存X密钥的文件中读数据
void readX(int row,const char *address){      //要读的行号(行号从1算起)，文件地址 
 
	ifstream ifile;
	ifile.open(address);
	char buff[3*n];
	int i;
	for(i=0;i<(row-1);i++){    //读并忽略前(row-1)行 
		ifile.getline(buff, 3*n, '\n'); 
	}
	for(i=0;i<(n-1);i++){      //读需要行的密钥 
		ifile.getline(buff, 50, ' ');  //读一个位
		X[i] = atoi(buff); 
	}
	ifile.getline(buff, 50, '\n'); 
	X[n-1] = atoi(buff);
	
	ifile.close();	
} 


//从存Y密钥的文件中读数据
void readY(int row,const char *address){      //要读的行号(行号从1算起)，文件地址 
 
	ifstream ifile;
	ifile.open(address);
	char buff[3*n];
	int i;
	for(i=0;i<(row-1);i++){    //读并忽略前(row-1)行 
		ifile.getline(buff, 3*n, '\n'); 
	}
	for(i=0;i<(n-1);i++){      //读需要行的密钥 
		ifile.getline(buff, 50, ' ');  //读一个位
		Y[i] = atoi(buff); 
	}
	ifile.getline(buff, 50, '\n'); 
	Y[n-1] = atoi(buff);
	
	ifile.close();	
} 

//返回n 
int getN(){
	return n;
}




