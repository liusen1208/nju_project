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

//边结构体 
struct Edge {
	int N;        //CN/VN 
	int numOfRL;  //Hr[][]（Hl[][]）该元素和 Hl[][numOfRL] (Hr[][numOfRL])相等   
	double q0,q1; //qij(0),qij(1);   //针对BP算法保留的项 
	double r0,r1; //rji(0),rji(1); 
};
vector<Edge> Hr1[m];         //H的行优先存储
vector<Edge> Hl1[n];         //H的列优先存储 
vector<Edge> Hr2[m];         
vector<Edge> Hl2[n];              

int Dv1[n] = { 0 };          //VN度分布 
int Dc1[m] = { 0 };          //CN度分布
int Dv2[n] = { 0 };          //VN度分布 
int Dc2[m] = { 0 };          //CN度分布

int Z1[m];                   //Z 
int Z2[m];                   //Z 

int X[n];                   //X
int Y[n];                   //Y 


/* 用来确定Hr和Hl中Edge结构体中的 numOfRL值 */    
void calNumOfRL(){                      
	int i,j,k,vn,cn,tmp;
	for(i=0;i<n;i++){      //遍历VN
	
		cn=Dv1[i];               //该VN度 (Dv已经赋初值)
		for(k=0;k<cn;k++){       //遍历与该VN相连的CN 
			j=Hl1[i][k].N;              //与VN i相连的CN j 
			tmp = 0;
			while(Hr1[j][tmp].N != i) tmp++;  //找出VN i为第几个与CN j相连VN 
			Hl1[i][k].numOfRL=tmp;      //求出numOfRL值           
		} 
		
		cn=Dv2[i];               //该VN度 (Dv已经赋初值)
		for(k=0;k<cn;k++){       //遍历与该VN相连的CN 
			j=Hl2[i][k].N;              //与VN i相连的CN j 
			tmp = 0;
			while(Hr2[j][tmp].N != i) tmp++;  //找出VN i为第几个与CN j相连VN 
			Hl2[i][k].numOfRL=tmp;      //求出numOfRL值           
		} 
	}
	
	for(j=0;j<m;j++){      //遍历CN
	
		vn=Dc1[j];                //该CN度 (Dc已经赋初值)
		for(k=0;k<vn;k++){       //遍历与该CN相连的VN 
			i=Hr1[j][k].N;              //与CN j相连的VN i
			tmp = 0;
			while(Hl1[i][tmp].N != j) tmp++;  //找出CN j为第几个与VN i相连CN 
			Hr1[j][k].numOfRL=tmp;      //求出numOfRL值           
		} 
		
		vn=Dc2[j];                //该CN度 (Dc已经赋初值)
		for(k=0;k<vn;k++){       //遍历与该CN相连的VN 
			i=Hr2[j][k].N;              //与CN j相连的VN i
			tmp = 0;
			while(Hl2[i][tmp].N != j) tmp++;  //找出CN j为第几个与VN i相连CN 
			Hr2[j][k].numOfRL=tmp;      //求出numOfRL值           
		} 
	}
}

/*从文本中读取H，并初始化 Dc 和 Dv*/
void readH(const char * address,int Num) {  //(文件地址，和读第几个矩阵)
	ifstream ifile;
	ifile.open(address);
	char buff[50];
	int row, column;
	Edge ed;
	
	cout<<"1.开始从 \""<<address<<"\" 读LDPC码！"<<endl;
	cout<<"    读取中..."<<endl; 
	
	if(Num==1){
		
		while (ifile.getline(buff, 50, ' ')) {//读出行 
			row = atoi(buff);                 	  //转成整型 
			ifile.getline(buff, 50, '\n');        //读出列 
			column = atoi(buff);                  //转成整型
			ed.N = column; Hr1[row].push_back(ed); //放到Hr,Hl 
			ed.N = row; Hl1[column].push_back(ed);
			Dc1[row]++;                            //该CN的度+1
			Dv1[column]++;                         //该VN的度+1
		}
	}else{
		while (ifile.getline(buff, 50, ' ')) {//读出行 
			row = atoi(buff);                 	  //转成整型 
			ifile.getline(buff, 50, '\n');        //读出列 
			column = atoi(buff);                  //转成整型
			ed.N = column; Hr2[row].push_back(ed); //放到Hr,Hl 
			ed.N = row; Hl2[column].push_back(ed);
			Dc2[row]++;                            //该CN的度+1
			Dv2[column]++;                         //该VN的度+1
		}
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

/*多矩阵BP算法，返回协商是否正常结束*/
int MBP2(double e,int iteraM){   //误码率，迭代次数 

    int i,j,k,kk,iteraN=0,tmpe; 
    double tmp0,tmp1,tmp3,tmp4;
    double tmp[11][2];
    bool isEqual;                 //判断H*Y是否等于Z 
    double P[n][2];               //VN的初始概率 
    double Q[n][2];               //VN的后验概率 
    int vn,cn;                    //暂存的VN/CN个数 
    
	//初始化Pi
	for(i=0;i<n;i++){       //对每个VN i 
		if(Y[i]==0){              //Yi取0 
			P[i][0]=1-e;               //Pi(0) 
			P[i][1]=e;	               //Pi(1)
		}else{                    //Yi取1 
			P[i][0]=e;                 //Pi(0)
			P[i][1]=1-e;	           //Pi(1)
		}
	} 
	
	//初始化qij 
	for(j=0;j<m;j++){       //遍历每个CN 
	
		vn=Dc1[j];                            //##################################
		for(i=0;i<vn;i++){//遍历与该CN相连VN 
		    tmpe=Hr1[j][i].N;                  //VN
		    Hr1[j][i].q0=P[tmpe][0];           //q tmpe_j(0)
			Hr1[j][i].q1=P[tmpe][1];           //q tmpe_j(1)
			k=Hr1[j][i].numOfRL;
			Hl1[tmpe][k].q0=Hr1[j][i].q0;       //q tmpe_j(0)
			Hl1[tmpe][k].q1=Hr1[j][i].q1;       //q tmpe_j(1)
		}
		
		vn=Dc2[j];                            //##################################
		for(i=0;i<vn;i++){//遍历与该CN相连VN 
		    tmpe=Hr2[j][i].N;                  //VN
		    Hr2[j][i].q0=P[tmpe][0];           //q tmpe_j(0)
			Hr2[j][i].q1=P[tmpe][1];           //q tmpe_j(1)
			k=Hr2[j][i].numOfRL;
			Hl2[tmpe][k].q0=Hr2[j][i].q0;       //q tmpe_j(0)
			Hl2[tmpe][k].q1=Hr2[j][i].q1;       //q tmpe_j(1)
		}
	}
	
	//处理CN/VN消息，译码判决，停止判断 
	do{
		iteraN++;                         //迭代次数加1 
		
		//处理CN消息 
		for(j=0;j<m;j++){                 //遍历每个CN
		
		    vn=Dc1[j];                           //与该CN相连VN数	
			kk=0;	
			for(tmp1=1.0,k=0;k<vn;k++)          //计算公式中一个乘积
				if(Hr1[j][k].q1!=0.5)  tmp1=tmp1*(1.0-2*Hr1[j][k].q1); 
				else kk++;		 
			for(i=0;i<vn;i++){                  //遍历这些VN 				
				tmpe=Hr1[j][i].N;                          //该VN 
				if(kk==0){
					tmp0=tmp1/(1.0-2*Hr1[j][i].q1); 
				}else{
					if(kk==1){
						if(Hr1[j][i].q1==0.5){
							tmp0=tmp1;
						}else{
							tmp0=0;
						}
					}else{
						tmp0=0;
					}
				}	
				Hr1[j][i].r0=(Z1[j]==1)?(0.5-0.5*tmp0):(0.5+0.5*tmp0);//r j_tmpe(0) 
				Hr1[j][i].r1=1.0-Hr1[j][i].r0;              //r j_tmpe(1) 
				k=Hr1[j][i].numOfRL;		                      
				Hl1[tmpe][k].r0=Hr1[j][i].r0;               //r j_tmpe(0)
				Hl1[tmpe][k].r1=Hr1[j][i].r1;               //r j_tmpe(1)		
			}
			
			vn=Dc2[j];                           //与该CN相连VN数
			kk=0;		
			for(tmp1=1.0,k=0;k<vn;k++)          //计算公式中一个乘积
				if(Hr2[j][k].q1!=0.5)  tmp1=tmp1*(1.0-2*Hr2[j][k].q1); 
				else kk++;			 
			for(i=0;i<vn;i++){                  //遍历这些VN 				
				tmpe=Hr2[j][i].N;                          //该VN 
				if(kk==0){
					tmp0=tmp1/(1.0-2*Hr2[j][i].q1); 
				}else{
					if(kk==1){
						if(Hr2[j][i].q1==0.5){
							tmp0=tmp1;
						}else{
							tmp0=0;
						}
					}else{
						tmp0=0;
					}
				}
				Hr2[j][i].r0=(Z2[j]==1)?(0.5-0.5*tmp0):(0.5+0.5*tmp0);//r j_tmpe(0) 
				Hr2[j][i].r1=1.0-Hr2[j][i].r0;              //r j_tmpe(1) 
				k=Hr2[j][i].numOfRL;		                      
				Hl2[tmpe][k].r0=Hr2[j][i].r0;               //r j_tmpe(0)
				Hl2[tmpe][k].r1=Hr2[j][i].r1;               //r j_tmpe(1)		
			}
		}
								
		//处理VN消息，译码判决 
		for(i=0;i<n;i++){                 //遍历每个VN
			
			//译码判决该VN 		
			cn=Dv1[i];                              //与该VN相连CN数  
			tmp[1][0]=1.0; tmp[1][1]=1.0;   //计算公式中两个乘积
			for(k=0;k<cn;k++){  
				tmp[1][0]=tmp[1][0]*Hl1[i][k].r0;
				tmp[1][1]=tmp[1][1]*Hl1[i][k].r1;
			}
			
			cn=Dv2[i];                     //与该VN相连CN数  
			tmp[2][0]=1.0; tmp[2][1]=1.0;   //计算公式中两个乘积
			for(k=0;k<cn;k++){  
				tmp[2][0]=tmp[2][0]*Hl2[i][k].r0;
				tmp[2][1]=tmp[2][1]*Hl2[i][k].r1;
			}	
			
			tmp0=P[i][0]*tmp[1][0]*tmp[2][0];
			tmp1=P[i][1]*tmp[1][1]*tmp[2][1];
			
			Q[i][0]=((tmp0+tmp1)==0.0)?0.5:((1/(tmp0+tmp1))*tmp0);//Qi(0)
			Q[i][1]=1.0-Q[i][0];                  //Qi(1)     //################################		
			if(Q[i][0]>0.5)  Y[i]=0;            //Qi(0)>0.5，Yi改为0 
			else if(Q[i][1]>0.5)  Y[i]=1;        //Qi(1)>0.5，Yi改为1			
						
			//处理该VN消息	
			cn=Dv1[i];
			for(j=0;j<cn;j++){
				tmp3=P[i][0]*tmp[1][0]/Hl1[i][j].r0;	tmp4=P[i][1]*tmp[1][1]/Hl1[i][j].r1;	
				Hl1[i][j].q0=((tmp3+tmp4)==0.0)?0.5:((1.0/(tmp3+tmp4))*tmp3);//q i_tmpe(0)	              
				Hl1[i][j].q1=1-Hl1[i][j].q0;            //q i_tmpe(1)
				tmpe=Hl1[i][j].N;  k=Hl1[i][j].numOfRL;                                 
				Hr1[tmpe][k].q0=Hl1[i][j].q0;           //q i_tmpe(0)
				Hr1[tmpe][k].q1=Hl1[i][j].q1;           //q i_tmpe(1)	
			}	
			
			cn=Dv2[i];
			for(j=0;j<cn;j++){
				tmp3=P[i][0]*tmp[2][0]/Hl2[i][j].r0;	tmp4=P[i][1]*tmp[2][1]/Hl2[i][j].r1;	
				Hl2[i][j].q0=((tmp3+tmp4)==0.0)?0.5:((1.0/(tmp3+tmp4))*tmp3);//q i_tmpe(0)	              
				Hl2[i][j].q1=1-Hl2[i][j].q0;            //q i_tmpe(1)
				tmpe=Hl2[i][j].N;  k=Hl2[i][j].numOfRL;                                 
				Hr2[tmpe][k].q0=Hl2[i][j].q0;           //q i_tmpe(0)
				Hr2[tmpe][k].q1=Hl2[i][j].q1;           //q i_tmpe(1)	
			}		                
		}	

		//停止判断
		isEqual = true;               //设H*Y=Z 
		for (j = 0;j<m;j++) {           //计算H*Y 
		
			vn=Dc1[j];		
			for (k=0,i = 0;i<vn;i++){       //遍历与 CN j相连的所有 VN 
				tmpe=Hr1[j][i].N;                  //对VN tmpe
				if(Y[tmpe]==1) k++;               //H第j行与Y对应位相乘为 1 的个数加1 
			}	
			k=((k%2)==1)?1:0;               //H第j行与Y相乘(为得到0/1结果，*/+分别用&&/^代替)
							
			if (k != Z1[j]) {                  //H*Y!=Z 
				isEqual = false;                 //标记H*Y!=Z 
				break;                         //停止剩余计算 
			}
		}
	} while(iteraN<iteraM&&!isEqual);//未达最大迭代，且H*Y!=Z 
	
	if(!isEqual){
		iteraN++;
	}
	return iteraN;       //返回迭代次数 
}

//获取Alice长度为m校验子Z
void getZ(){
	int i,j,k,vn;
	
	for (j = 0;j<m;j++) {           //计算H*X		
	
		vn=Dc1[j]; k=0;
		for(i=0;i<vn;i++){          //遍历与 CN j相连的所有 VN 
			if(X[Hr1[j][i].N]==1) k++;     //H第j行与Y对应位相乘为 1 的个数加1 
		}	
		Z1[j]=((k%2)==1)?1:0;        //H第j行与Y相乘(为得到0/1结果，*/+分别用&&/^代替)
		
		vn=Dc2[j]; k=0;
		for(i=0;i<vn;i++){          //遍历与 CN j相连的所有 VN 
			if(X[Hr2[j][i].N]==1) k++;     //H第j行与Y对应位相乘为 1 的个数加1 
		}	
		Z2[j]=((k%2)==1)?1:0;        //H第j行与Y相乘(为得到0/1结果，*/+分别用&&/^代替)
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
	for(i=0;i<errorN;i++){  //调控错误数    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
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

int getN(){
	return n;
}






