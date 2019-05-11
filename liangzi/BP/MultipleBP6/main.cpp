#include "Reconciliation.h" 
#include <iostream>
#include <vector>
#include <fstream>
#include <ctime> 
#include <cstdlib>
#include <cmath>
using namespace std;

//ע����&&&&&&&&&&&&&&&&&&��Ϊ������Ҫ�Ķ��Ĳ���λ�� 

int main(){	

	clock_t a,b;
	double c=0.0,d;
	double errorRate=0.0246;                                //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
	int afterWrong=0; 
  	int iteratBound=100;            //LLRBP�㷨�ĵ�������  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	int nn=100;                     //ѭ������             //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	int iteratN=0;                 //�����ܴ��� 
	int successN=0;                //Э�̳ɹ����� 
	int tmp,tmp1;	     
	
	cout<<"��ӭʹ�� MBP Э�̲��Գ���"<<endl<<endl; 
	 	
	readH("1.txt",1);                                //���ı��ж�ȡУ�����H������ʼ�� Dc �� Dv
	readH("2.txt",2);
	readH("3.txt",3);
	readH("4.txt",4);
	readH("5.txt",5);
	readH("6.txt",6);
	
	cout<<"2.LDPC���м���������..."<<endl;
	calNumOfRL();                                   //ȷ��Hr��Hl��Edge�ṹ���е� numOfRLֵ  
	cout<<"  LDPC���м���������ϣ�"<<endl<<endl;
	srand((unsigned)time(NULL));                    //������������� 
	
//	for(int i=0;i<nn;i++){          
//		cout<<"��"<<i+1<<"����ԿЭ�̲���......"<<endl; 
//		cout<<"3.Alice��ʼ�Ʊ��뷢��������Կ..."<<endl;
//		cout<<"    ������..."<<endl;
//		randomX();                   //���ģ��X 
//		cout<<"  Alice��������ϣ�"<<endl;
//		cout<<"  Bob����ɽ��գ�"<<endl;  		
//		cout<<"4.��ɶԻ���"<<endl;
//		cout<<"5.��ɲ������ƣ�"<<endl; 
//		
//		cout<<"6.Alice����У���ӣ�"<<endl; 
//		getZ();                      //��ȡAlice����ΪmУ����Z   
//		randomY(errorRate);                   //��Xģ���Y  
//		
//		cout<<"7.Bob��ʼ��ԿЭ��..."<<endl; 
//		a=clock();	
//		tmp=MBP6(errorRate,iteratBound);    //BP�㷨 
//		b=clock();
//		
//		if(tmp<=iteratBound){   
//		   
//			cout<<"    ����Э�̳ɹ���"<<endl; 
//			successN++;   
//		}else{
//			
//			cout<<"    ����Э��ʧ�ܣ�"<<endl; 
//			tmp--;
//		} 
//		iteratN+=tmp;
//		cout<<"    ��"<<tmp<<"�ε�������"<<endl;
//		
//		tmp1=Test();
//		if(tmp1!=0){
//			cout<<"    ʣ���������"<<tmp1<<endl;
//		}
//		afterWrong+=tmp1;
//		
//		d=(double)(b-a)/CLOCKS_PER_SEC;
//		c=c+d;
//		cout<<"    ʱ�䣺"<<d<<"s"<<endl<<endl; 	
//	}  
	
	
	
	//createXY(nn,errorRate);      //����nn��X��Y����д���ļ� 
	for(int i=0;i<nn;i++){          
		cout<<"��"<<i+1<<"����ԿЭ�̲���......"<<endl; 
		cout<<"3.Alice��ʼ�Ʊ��뷢��������Կ..."<<endl;
		cout<<"    ������..."<<endl;
		
		readX(i+1,"X_100_0.0246.txt");                   //��X           //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		
		cout<<"  Alice��������ϣ�"<<endl;
		cout<<"  Bob����ɽ��գ�"<<endl;  		
		cout<<"4.��ɶԻ���"<<endl;
		cout<<"5.��ɲ������ƣ�"<<endl; 
		
		cout<<"6.Alice����У���ӣ�"<<endl; 
		getZ();                      //��ȡAlice����ΪmУ����Z   
		readY(i+1,"Y_100_0.0246.txt");                   //��Y           //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		
		cout<<"7.Bob��ʼ��ԿЭ��..."<<endl; 
		a=clock();
		tmp=MBP6(errorRate,iteratBound);    //BP�㷨 	 
		b=clock();
		
		if(tmp<=iteratBound){   
		   
			cout<<"    ����Э�̳ɹ���"<<endl; 
			successN++;   
		}else{
			
			cout<<"    ����Э��ʧ�ܣ�"<<endl; 
			tmp--;
		} 
		iteratN+=tmp;
		cout<<"    ��"<<tmp<<"�ε�������"<<endl;
		
		tmp1=Test();
		if(tmp1!=0){
			cout<<"    ʣ���������"<<tmp1<<endl;
		}
		afterWrong+=tmp1;
		
		d=(double)(b-a)/CLOCKS_PER_SEC;
		c=c+d;
		cout<<"    ʱ�䣺"<<d<<"s"<<endl<<endl; 	
	}   		 
	
	cout<<"���Դ�����"<<nn<<endl;
	cout<<"�����ʣ�"<<errorRate<<endl; 
	cout<<"ƽ��ʱ�䣺"<<c/nn<<endl;	
	cout<<"ƽ������������"<<iteratN*1.0/nn<<endl;
	cout<<"��������ʣ�"<<afterWrong*1.0/(getN()*nn)<<endl;      
	cout<<"Э�̳ɹ�������"<<successN<<endl;
	cout<<"Э�̳ɹ��ʣ�"<<successN*1.0/nn<<endl; 

}













