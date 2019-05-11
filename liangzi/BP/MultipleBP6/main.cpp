#include "Reconciliation.h" 
#include <iostream>
#include <vector>
#include <fstream>
#include <ctime> 
#include <cstdlib>
#include <cmath>
using namespace std;

//注：“&&&&&&&&&&&&&&&&&&”为可能需要改动的参数位置 

int main(){	

	clock_t a,b;
	double c=0.0,d;
	double errorRate=0.0246;                                //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
	int afterWrong=0; 
  	int iteratBound=100;            //LLRBP算法的迭代上限  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	int nn=100;                     //循环次数             //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	int iteratN=0;                 //迭代总次数 
	int successN=0;                //协商成功次数 
	int tmp,tmp1;	     
	
	cout<<"欢迎使用 MBP 协商测试程序！"<<endl<<endl; 
	 	
	readH("1.txt",1);                                //从文本中读取校验矩阵H，并初始化 Dc 和 Dv
	readH("2.txt",2);
	readH("3.txt",3);
	readH("4.txt",4);
	readH("5.txt",5);
	readH("6.txt",6);
	
	cout<<"2.LDPC码中间量计算中..."<<endl;
	calNumOfRL();                                   //确定Hr和Hl中Edge结构体中的 numOfRL值  
	cout<<"  LDPC码中间量计算完毕！"<<endl<<endl;
	srand((unsigned)time(NULL));                    //设置随机数种子 
	
//	for(int i=0;i<nn;i++){          
//		cout<<"第"<<i+1<<"次密钥协商测试......"<<endl; 
//		cout<<"3.Alice开始制备与发送量子密钥..."<<endl;
//		cout<<"    发送中..."<<endl;
//		randomX();                   //随机模拟X 
//		cout<<"  Alice方发送完毕！"<<endl;
//		cout<<"  Bob方完成接收！"<<endl;  		
//		cout<<"4.完成对基！"<<endl;
//		cout<<"5.完成参数估计！"<<endl; 
//		
//		cout<<"6.Alice发送校验子！"<<endl; 
//		getZ();                      //获取Alice长度为m校验子Z   
//		randomY(errorRate);                   //由X模拟出Y  
//		
//		cout<<"7.Bob开始密钥协商..."<<endl; 
//		a=clock();	
//		tmp=MBP6(errorRate,iteratBound);    //BP算法 
//		b=clock();
//		
//		if(tmp<=iteratBound){   
//		   
//			cout<<"    本次协商成功！"<<endl; 
//			successN++;   
//		}else{
//			
//			cout<<"    本次协商失败！"<<endl; 
//			tmp--;
//		} 
//		iteratN+=tmp;
//		cout<<"    第"<<tmp<<"次迭代结束"<<endl;
//		
//		tmp1=Test();
//		if(tmp1!=0){
//			cout<<"    剩余错误数："<<tmp1<<endl;
//		}
//		afterWrong+=tmp1;
//		
//		d=(double)(b-a)/CLOCKS_PER_SEC;
//		c=c+d;
//		cout<<"    时间："<<d<<"s"<<endl<<endl; 	
//	}  
	
	
	
	//createXY(nn,errorRate);      //创造nn条X和Y，并写入文件 
	for(int i=0;i<nn;i++){          
		cout<<"第"<<i+1<<"次密钥协商测试......"<<endl; 
		cout<<"3.Alice开始制备与发送量子密钥..."<<endl;
		cout<<"    发送中..."<<endl;
		
		readX(i+1,"X_100_0.0246.txt");                   //读X           //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		
		cout<<"  Alice方发送完毕！"<<endl;
		cout<<"  Bob方完成接收！"<<endl;  		
		cout<<"4.完成对基！"<<endl;
		cout<<"5.完成参数估计！"<<endl; 
		
		cout<<"6.Alice发送校验子！"<<endl; 
		getZ();                      //获取Alice长度为m校验子Z   
		readY(i+1,"Y_100_0.0246.txt");                   //读Y           //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		
		cout<<"7.Bob开始密钥协商..."<<endl; 
		a=clock();
		tmp=MBP6(errorRate,iteratBound);    //BP算法 	 
		b=clock();
		
		if(tmp<=iteratBound){   
		   
			cout<<"    本次协商成功！"<<endl; 
			successN++;   
		}else{
			
			cout<<"    本次协商失败！"<<endl; 
			tmp--;
		} 
		iteratN+=tmp;
		cout<<"    第"<<tmp<<"次迭代结束"<<endl;
		
		tmp1=Test();
		if(tmp1!=0){
			cout<<"    剩余错误数："<<tmp1<<endl;
		}
		afterWrong+=tmp1;
		
		d=(double)(b-a)/CLOCKS_PER_SEC;
		c=c+d;
		cout<<"    时间："<<d<<"s"<<endl<<endl; 	
	}   		 
	
	cout<<"测试次数："<<nn<<endl;
	cout<<"错误率："<<errorRate<<endl; 
	cout<<"平均时间："<<c/nn<<endl;	
	cout<<"平均迭代次数："<<iteratN*1.0/nn<<endl;
	cout<<"纠后错误率："<<afterWrong*1.0/(getN()*nn)<<endl;      
	cout<<"协商成功次数："<<successN<<endl;
	cout<<"协商成功率："<<successN*1.0/nn<<endl; 

}













