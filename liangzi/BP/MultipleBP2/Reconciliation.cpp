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
#define n 10000                //H�Ĺ�ģ��m*n       //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
#define m 2000                                      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
using namespace std;

/*
�����ڵ�(variance node),У��ڵ㣨check node���ֱ���ΪVN,CN ��
У�������H��ʾ��   Alice/Bobɸ����Կ ��X/Y��    AliceУ���ӣ�Z
*/

//�߽ṹ�� 
struct Edge {
	int N;        //CN/VN 
	int numOfRL;  //Hr[][]��Hl[][]����Ԫ�غ� Hl[][numOfRL] (Hr[][numOfRL])���   
	double q0,q1; //qij(0),qij(1);   //���BP�㷨�������� 
	double r0,r1; //rji(0),rji(1); 
};
vector<Edge> Hr1[m];         //H�������ȴ洢
vector<Edge> Hl1[n];         //H�������ȴ洢 
vector<Edge> Hr2[m];         
vector<Edge> Hl2[n];              

int Dv1[n] = { 0 };          //VN�ȷֲ� 
int Dc1[m] = { 0 };          //CN�ȷֲ�
int Dv2[n] = { 0 };          //VN�ȷֲ� 
int Dc2[m] = { 0 };          //CN�ȷֲ�

int Z1[m];                   //Z 
int Z2[m];                   //Z 

int X[n];                   //X
int Y[n];                   //Y 


/* ����ȷ��Hr��Hl��Edge�ṹ���е� numOfRLֵ */    
void calNumOfRL(){                      
	int i,j,k,vn,cn,tmp;
	for(i=0;i<n;i++){      //����VN
	
		cn=Dv1[i];               //��VN�� (Dv�Ѿ�����ֵ)
		for(k=0;k<cn;k++){       //�������VN������CN 
			j=Hl1[i][k].N;              //��VN i������CN j 
			tmp = 0;
			while(Hr1[j][tmp].N != i) tmp++;  //�ҳ�VN iΪ�ڼ�����CN j����VN 
			Hl1[i][k].numOfRL=tmp;      //���numOfRLֵ           
		} 
		
		cn=Dv2[i];               //��VN�� (Dv�Ѿ�����ֵ)
		for(k=0;k<cn;k++){       //�������VN������CN 
			j=Hl2[i][k].N;              //��VN i������CN j 
			tmp = 0;
			while(Hr2[j][tmp].N != i) tmp++;  //�ҳ�VN iΪ�ڼ�����CN j����VN 
			Hl2[i][k].numOfRL=tmp;      //���numOfRLֵ           
		} 
	}
	
	for(j=0;j<m;j++){      //����CN
	
		vn=Dc1[j];                //��CN�� (Dc�Ѿ�����ֵ)
		for(k=0;k<vn;k++){       //�������CN������VN 
			i=Hr1[j][k].N;              //��CN j������VN i
			tmp = 0;
			while(Hl1[i][tmp].N != j) tmp++;  //�ҳ�CN jΪ�ڼ�����VN i����CN 
			Hr1[j][k].numOfRL=tmp;      //���numOfRLֵ           
		} 
		
		vn=Dc2[j];                //��CN�� (Dc�Ѿ�����ֵ)
		for(k=0;k<vn;k++){       //�������CN������VN 
			i=Hr2[j][k].N;              //��CN j������VN i
			tmp = 0;
			while(Hl2[i][tmp].N != j) tmp++;  //�ҳ�CN jΪ�ڼ�����VN i����CN 
			Hr2[j][k].numOfRL=tmp;      //���numOfRLֵ           
		} 
	}
}

/*���ı��ж�ȡH������ʼ�� Dc �� Dv*/
void readH(const char * address,int Num) {  //(�ļ���ַ���Ͷ��ڼ�������)
	ifstream ifile;
	ifile.open(address);
	char buff[50];
	int row, column;
	Edge ed;
	
	cout<<"1.��ʼ�� \""<<address<<"\" ��LDPC�룡"<<endl;
	cout<<"    ��ȡ��..."<<endl; 
	
	if(Num==1){
		
		while (ifile.getline(buff, 50, ' ')) {//������ 
			row = atoi(buff);                 	  //ת������ 
			ifile.getline(buff, 50, '\n');        //������ 
			column = atoi(buff);                  //ת������
			ed.N = column; Hr1[row].push_back(ed); //�ŵ�Hr,Hl 
			ed.N = row; Hl1[column].push_back(ed);
			Dc1[row]++;                            //��CN�Ķ�+1
			Dv1[column]++;                         //��VN�Ķ�+1
		}
	}else{
		while (ifile.getline(buff, 50, ' ')) {//������ 
			row = atoi(buff);                 	  //ת������ 
			ifile.getline(buff, 50, '\n');        //������ 
			column = atoi(buff);                  //ת������
			ed.N = column; Hr2[row].push_back(ed); //�ŵ�Hr,Hl 
			ed.N = row; Hl2[column].push_back(ed);
			Dc2[row]++;                            //��CN�Ķ�+1
			Dv2[column]++;                         //��VN�Ķ�+1
		}
	}
	
	ifile.close();	
	cout<<"  LDPC���ȡ������"<<endl<<endl;
}

//����X��Y�ж��ٲ�ͬ 
int Test() {  //��һ����������Կ������
	int sum = 0, i;
	for (i = 0;i < n;i++) {
		if (X[i] != Y[i])
			sum++;
	}
	return sum;
}

/*�����BP�㷨������Э���Ƿ���������*/
int MBP2(double e,int iteraM){   //�����ʣ��������� 

    int i,j,k,kk,iteraN=0,tmpe; 
    double tmp0,tmp1,tmp3,tmp4;
    double tmp[11][2];
    bool isEqual;                 //�ж�H*Y�Ƿ����Z 
    double P[n][2];               //VN�ĳ�ʼ���� 
    double Q[n][2];               //VN�ĺ������ 
    int vn,cn;                    //�ݴ��VN/CN���� 
    
	//��ʼ��Pi
	for(i=0;i<n;i++){       //��ÿ��VN i 
		if(Y[i]==0){              //Yiȡ0 
			P[i][0]=1-e;               //Pi(0) 
			P[i][1]=e;	               //Pi(1)
		}else{                    //Yiȡ1 
			P[i][0]=e;                 //Pi(0)
			P[i][1]=1-e;	           //Pi(1)
		}
	} 
	
	//��ʼ��qij 
	for(j=0;j<m;j++){       //����ÿ��CN 
	
		vn=Dc1[j];                            //##################################
		for(i=0;i<vn;i++){//�������CN����VN 
		    tmpe=Hr1[j][i].N;                  //VN
		    Hr1[j][i].q0=P[tmpe][0];           //q tmpe_j(0)
			Hr1[j][i].q1=P[tmpe][1];           //q tmpe_j(1)
			k=Hr1[j][i].numOfRL;
			Hl1[tmpe][k].q0=Hr1[j][i].q0;       //q tmpe_j(0)
			Hl1[tmpe][k].q1=Hr1[j][i].q1;       //q tmpe_j(1)
		}
		
		vn=Dc2[j];                            //##################################
		for(i=0;i<vn;i++){//�������CN����VN 
		    tmpe=Hr2[j][i].N;                  //VN
		    Hr2[j][i].q0=P[tmpe][0];           //q tmpe_j(0)
			Hr2[j][i].q1=P[tmpe][1];           //q tmpe_j(1)
			k=Hr2[j][i].numOfRL;
			Hl2[tmpe][k].q0=Hr2[j][i].q0;       //q tmpe_j(0)
			Hl2[tmpe][k].q1=Hr2[j][i].q1;       //q tmpe_j(1)
		}
	}
	
	//����CN/VN��Ϣ�������о���ֹͣ�ж� 
	do{
		iteraN++;                         //����������1 
		
		//����CN��Ϣ 
		for(j=0;j<m;j++){                 //����ÿ��CN
		
		    vn=Dc1[j];                           //���CN����VN��	
			kk=0;	
			for(tmp1=1.0,k=0;k<vn;k++)          //���㹫ʽ��һ���˻�
				if(Hr1[j][k].q1!=0.5)  tmp1=tmp1*(1.0-2*Hr1[j][k].q1); 
				else kk++;		 
			for(i=0;i<vn;i++){                  //������ЩVN 				
				tmpe=Hr1[j][i].N;                          //��VN 
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
			
			vn=Dc2[j];                           //���CN����VN��
			kk=0;		
			for(tmp1=1.0,k=0;k<vn;k++)          //���㹫ʽ��һ���˻�
				if(Hr2[j][k].q1!=0.5)  tmp1=tmp1*(1.0-2*Hr2[j][k].q1); 
				else kk++;			 
			for(i=0;i<vn;i++){                  //������ЩVN 				
				tmpe=Hr2[j][i].N;                          //��VN 
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
								
		//����VN��Ϣ�������о� 
		for(i=0;i<n;i++){                 //����ÿ��VN
			
			//�����о���VN 		
			cn=Dv1[i];                              //���VN����CN��  
			tmp[1][0]=1.0; tmp[1][1]=1.0;   //���㹫ʽ�������˻�
			for(k=0;k<cn;k++){  
				tmp[1][0]=tmp[1][0]*Hl1[i][k].r0;
				tmp[1][1]=tmp[1][1]*Hl1[i][k].r1;
			}
			
			cn=Dv2[i];                     //���VN����CN��  
			tmp[2][0]=1.0; tmp[2][1]=1.0;   //���㹫ʽ�������˻�
			for(k=0;k<cn;k++){  
				tmp[2][0]=tmp[2][0]*Hl2[i][k].r0;
				tmp[2][1]=tmp[2][1]*Hl2[i][k].r1;
			}	
			
			tmp0=P[i][0]*tmp[1][0]*tmp[2][0];
			tmp1=P[i][1]*tmp[1][1]*tmp[2][1];
			
			Q[i][0]=((tmp0+tmp1)==0.0)?0.5:((1/(tmp0+tmp1))*tmp0);//Qi(0)
			Q[i][1]=1.0-Q[i][0];                  //Qi(1)     //################################		
			if(Q[i][0]>0.5)  Y[i]=0;            //Qi(0)>0.5��Yi��Ϊ0 
			else if(Q[i][1]>0.5)  Y[i]=1;        //Qi(1)>0.5��Yi��Ϊ1			
						
			//�����VN��Ϣ	
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

		//ֹͣ�ж�
		isEqual = true;               //��H*Y=Z 
		for (j = 0;j<m;j++) {           //����H*Y 
		
			vn=Dc1[j];		
			for (k=0,i = 0;i<vn;i++){       //������ CN j���������� VN 
				tmpe=Hr1[j][i].N;                  //��VN tmpe
				if(Y[tmpe]==1) k++;               //H��j����Y��Ӧλ���Ϊ 1 �ĸ�����1 
			}	
			k=((k%2)==1)?1:0;               //H��j����Y���(Ϊ�õ�0/1�����*/+�ֱ���&&/^����)
							
			if (k != Z1[j]) {                  //H*Y!=Z 
				isEqual = false;                 //���H*Y!=Z 
				break;                         //ֹͣʣ����� 
			}
		}
	} while(iteraN<iteraM&&!isEqual);//δ������������H*Y!=Z 
	
	if(!isEqual){
		iteraN++;
	}
	return iteraN;       //���ص������� 
}

//��ȡAlice����ΪmУ����Z
void getZ(){
	int i,j,k,vn;
	
	for (j = 0;j<m;j++) {           //����H*X		
	
		vn=Dc1[j]; k=0;
		for(i=0;i<vn;i++){          //������ CN j���������� VN 
			if(X[Hr1[j][i].N]==1) k++;     //H��j����Y��Ӧλ���Ϊ 1 �ĸ�����1 
		}	
		Z1[j]=((k%2)==1)?1:0;        //H��j����Y���(Ϊ�õ�0/1�����*/+�ֱ���&&/^����)
		
		vn=Dc2[j]; k=0;
		for(i=0;i<vn;i++){          //������ CN j���������� VN 
			if(X[Hr2[j][i].N]==1) k++;     //H��j����Y��Ӧλ���Ϊ 1 �ĸ�����1 
		}	
		Z2[j]=((k%2)==1)?1:0;        //H��j����Y���(Ϊ�õ�0/1�����*/+�ֱ���&&/^����)
	}	
} 

//ģ��X
void randomX(){
	int i;
	for(i=0;i<n;i++){
		X[i]=rand() % 2;
	}
} 

//��Xģ���Y 
void randomY(double e){
	int errorN=round(n*e);
	int i,k;
	int mark[n];
	for(i=0;i<n;i++){
		Y[i]=X[i];
		mark[i]=0;
	}
	for(i=0;i<errorN;i++){  //���ش�����    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
		k=rand()%n;
		if(mark[k]==1){
			i--;
			continue;
		}
		Y[k]=(Y[k]==0)?1:0;
		mark[k]=1;
	}
}

//�����ʺʹ����ʣ�����Э��Ч��f
double calF(double e){
	double H=-e*(log(e)/log(2))-(1-e)*(log(1-e)/log(2));
	return m*1.0/n/H;
} 


//����ָ��������ָ�������ʵ���ԿX��Y����д�������ļ���
void createXY(int num,double e){
	int k,kk,i;
	int errorN=round(n*e);
	int mark[n];
	ofstream ofile1,ofile2;
	ofile1.open("X.txt");
	ofile2.open("Y.txt");
	
	for(kk=0;kk<num;kk++){               //����num����Կ 
		for(i=0;i<(n-1);i++){               //������Կ��ǰ(n-1)λ 
			X[i]=rand() % 2;                      //����X��д���ļ� 
			ofile1<<X[i]<<" ";
		}
		X[n-1]=rand() % 2;                  //������Կ�����һλ            
		ofile1<<X[n-1]<<endl;
		
		for(i=0;i<n;i++){
			Y[i]=X[i];
			mark[i]=0;
		}
		for(i=0;i<errorN;i++){          //���ش�����    
			k=rand()%n;
			if(mark[k]==1){
				i--;
				continue;
			}
			Y[k]=(Y[k]==0)?1:0;
			mark[k]=1;
		}
		for(i=0;i<(n-1);i++){        //Yд���ļ� 
			ofile2<<Y[i]<<" ";
		} 
		ofile2<<Y[n-1]<<endl;
		
	}
	
	ofile1.close();
	ofile2.close();
} 


//�Ӵ�X��Կ���ļ��ж�����
void readX(int row,const char *address){      //Ҫ�����к�(�кŴ�1����)���ļ���ַ 
 
	ifstream ifile;
	ifile.open(address);
	char buff[3*n];
	int i;
	for(i=0;i<(row-1);i++){    //��������ǰ(row-1)�� 
		ifile.getline(buff, 3*n, '\n'); 
	}
	for(i=0;i<(n-1);i++){      //����Ҫ�е���Կ 
		ifile.getline(buff, 50, ' ');  //��һ��λ
		X[i] = atoi(buff); 
	}
	ifile.getline(buff, 50, '\n'); 
	X[n-1] = atoi(buff);
	
	ifile.close();	
} 


//�Ӵ�Y��Կ���ļ��ж�����
void readY(int row,const char *address){      //Ҫ�����к�(�кŴ�1����)���ļ���ַ 
 
	ifstream ifile;
	ifile.open(address);
	char buff[3*n];
	int i;
	for(i=0;i<(row-1);i++){    //��������ǰ(row-1)�� 
		ifile.getline(buff, 3*n, '\n'); 
	}
	for(i=0;i<(n-1);i++){      //����Ҫ�е���Կ 
		ifile.getline(buff, 50, ' ');  //��һ��λ
		Y[i] = atoi(buff); 
	}
	ifile.getline(buff, 50, '\n'); 
	Y[n-1] = atoi(buff);
	
	ifile.close();	
} 

int getN(){
	return n;
}






