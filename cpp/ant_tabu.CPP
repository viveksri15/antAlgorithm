#include<time.h>
#include<math.h>
#include<iostream.h>
#include<conio.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<fstream.h>

class ants
{
	int k,m,n,q;
	long double nmax,t,v,y,nmin;
	long double **km;
	int **pm,*btour;
	//km=phermone maintainance
	//pm=part machine matrix
	int *generate();
	float evaluate(int *);
	void updatepher(int *,float);
	void evaporate();
	void show_pm(int *);
	public:
	ants();
	~ants()
	{
		delete km,pm;
	}
	void start();
};

int *ants::generate()
{
	int *tour=new int[m+n];
	for(int i=0;i<m+n;i++)
	{
		long double sum=0;
		/**********normallizing*************/
		long double *pher=new long double[k];
		for(int j=0;j<k;j++)
		sum=sum+km[i][j];
		if(sum==0)continue;
		for(j=0;j<k;j++)
			pher[j]=km[i][j]/sum;
		int max=0;
		/********choosing apparant sol********/
		for(j=1;j<k;j++)
		{
			if(pher[max]<pher[j])
				max=j;
			if(pher[max]==pher[j])
			{
			float ran1=random(4);
			if(ran1>=2)max=j;
			}
		}
		/********choosing random sol**********/
		int ranch=0;
		for(j=1;j<k;j++)
			pher[j]+=pher[j-1];
		float ran1=(float)(random(pher[k-1]*1000))/1000;
		for(j=1;j<k;j++)
		{
			if(ran1>=0 && ran1<pher[0])
			{
				ranch=0;
				break;
			}
			if(ran1>=pher[j-1] && ran1<pher[j])
			{
				ranch=j;
				break;
			}
		}
		/***************************************/
		float ran=(float)(random(1000))/1000;
		if(ran>v)
			max=ranch;
		tour[i]=max;
		delete pher;
	}
	float efft1=evaluate(tour);
	int **tabu=new int*[k];
	for(int u2=0;u2<k;u2++)
		tabu[u2]=new int[2];
	int fill=0;
	int *toure=new int[m+n];
	for(u2=0;u2<m+n;u2++)
		toure[u2]=tour[u2];
	for(int u=0;u<k+1;u++)
	{
		int ran1=random(m+n);
		int ran2=random(k);
		int flag=0;
		for(int u1=0;u1<fill;u1++)
		{
			if(tabu[u1][0]==ran1 && tabu[u1][1]==ran2)
			{
				flag=1;
				break;
			}
		}
		if(flag==1){u--;continue;}
		tabu[fill][0]=ran1;
		tabu[fill][1]=ran2;
		fill++;
		int *tour1=new int[m+n];
		for(u1=0;u1<m+n;u1++)
			tour1[u1]=tour[u1];
		tour1[ran1]=ran2;
		if(efft1<evaluate(tour1))
			{for(u1=0;u1<m+n;u1++)
				toure[u1]=tour1[u1];
			efft1=evaluate(toure);
			}
		delete tour1;
	}
	delete tour,tabu;
	return toure;
}

float ants::evaluate(int *tour)
{
	/*** calculate diagonallized entries****/
	int q1=0,n1=0,n2=0;
	for(int i=m;i<m+n;i++)
	{
		for(int j=0;j<m;j++)
		{
			if(tour[i]==tour[j])
			{q1++;
			if(pm[j][i-m]==1)n1++;
			}
			else if(pm[j][i-m]==0)n2++;
		}
	}
	if(q1==0 || m*n==q1)return 0;
	long double qq=(long double)q1/(m*n);
	long double nn1=(long double)n1/q1,nn2=(long double)n2/(m*n-q1);
	long double eff=qq*nn1+(1-qq)*nn2;
	return eff;
}

void ants::updatepher(int *tour,float eff)
{
	for(int i=0;i<m+n;i++)
	{
		int tt=tour[i];
		km[i][tt]+=eff;
	}
}

void ants::evaporate()
{
	for(int i=0;i<m+n;i++)
		for(int j=0;j<k;j++)
			km[i][j]=km[i][j]*(1-y);
}

void ants::start()
{
	clrscr();
	randomize();
	fstream fin("tabu.txt",ios::out);
	nmin=t;nmax=0;
	for(int i=0;i<q;i++)
	{
		int *tour=generate();
		float eff1=evaluate(tour);
		//cout<<eff1<<' '<<i<<' ';
		//getch();
		fin<<eff1<<'\t';
		for(int xi=0;xi<m+n;xi++)
			fin<<btour[xi]<<' ';
		fin<<endl;
		if(eff1>nmax && eff1>nmin)
		{
			//cout<<eff1<<endl;
			//fin<<eff1<<endl;
			nmax=eff1;
			for(int xi=0;xi<m+n;xi++)
				btour[xi]=tour[xi];
			updatepher(tour,eff1);
		}
		evaporate();
		//cout<<endl;
		delete tour;
	}
	cout<<"best tour"<<endl;
	for(int xi=0;xi<m+n;xi++)
		cout<<btour[xi]<<' ';
	cout<<"max eff= : "<<nmax<<endl;
	show_pm(btour);
	fin<<nmax;
	fin.close();
}

void ants::show_pm(int *btour)
{
	int *mac=new int[m];
	int *par=new int[n];
	int pos=0;
	for(int i=0;i<k;i++)
	{
		for(int j=0;j<m;j++)
			if(btour[j]==i)
				mac[pos++]=j;
	}
	pos=0;
	for(i=0;i<k;i++)
	{
		for(int j=m;j<m+n;j++)
			if(btour[j]==i)
				par[pos++]=j-m;
	}
/*	cout<<endl;
	for(i=0;i<m;i++)
		cout<<mac[i]<<' ';
	cout<<endl;
	for(i=0;i<n;i++)
		cout<<par[i]<<' ';*/
	cout<<"generated machine-part matrix : ";
	cout<<endl;
	for(i=0;i<m;i++)
		{
		for(int j=0;j<n;j++)
			cout<<pm[mac[i]][par[j]]<<' ';
		cout<<endl;
		}

}

ants::ants()
{
	clrscr();
	fstream fin("matrix.txt",ios::in);
	fin>>k;
	fin>>m;
	fin>>n;
	fin>>t;
	fin>>v;
	fin>>y;
	fin>>q;
	pm=new int*[m];
	for(int i1=0;i1<m;i1++)
		pm[i1]=new int[n];
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			fin>>pm[i][j];
	btour=new int[m+n];
/*	for(i=0;i<m;i++)
	{
	for(int j=0;j<n;j++)
		cout<<pm[i][j]<<" ";
	cout<<endl;
	}*/
	fin.close();
	km=new long double*[m+n];
	for(i1=0;i1<m+n;i1++)
		km[i1]=new long double[k];
	for(i1=0;i1<m+n;i1++)
		btour[i1]=0;
	for(i=0;i<m+n;i++)
		for(j=0;j<k;j++)
			km[i][j]=1;
	//reading from files for input complete
}

void main()
{
	ants ant;
	ant.start();
	getch();
}
