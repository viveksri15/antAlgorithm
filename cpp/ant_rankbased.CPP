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
	int k,m,n,q,an,w;
	long double nmax,t,v,y,nmin;
	long double **km,**km1;
	int **pm,**rt,*btour;
	//km=phermone maintainance
	//pm=part machine matrix
	int *generate();
	float evaluate(int *);
	void updatepher(int *,float,int);
	void insert(int *,float);
	void show_pm(int *);
	void evaporate();
	public:
	ants();
	~ants()
	{
		delete km,pm,km1;
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
	return tour;
}

float ants::evaluate(int *tour)
{
	/*** calculate diagonallized entries****/
	if(tour[0]==-1)return -1;
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
	if(q1==0 || m*n-q1==0)return 0;
	long double qq=(long double)q1/(m*n);
	long double nn1=(long double)n1/q1,nn2=(long double)n2/(m*n-q1);
	long double eff=qq*nn1+(1-qq)*nn2;
	return eff;
}

void ants::updatepher(int *tour,float eff,int j)
{
	for(int i=0;i<m+n;i++)
	{
		int tt=tour[i];
		km1[i][tt]+=(w-j)*eff;
	}
}

void ants::evaporate()
{
	for(int i=0;i<m+n;i++)
		for(int j=0;j<k;j++)
			km1[i][j]*=(1-y);
}

void ants::start()
{
	clrscr();
	randomize();
	fstream fin("rank.txt",ios::out);
	nmax=nmin=t;
	for(int i=0;i<q/an;i++)
	{
		for(int j=0;j<m+n;j++)
			for(int j1=0;j1<k;j1++)
				km[j][j1]=km1[j][j1];
		for(j=0;j<an;j++)
		{
			int *tour=generate();
			float eff1=evaluate(tour);
			//cout<<eff1<<endl;
			fin<<eff1<<'\t';
			for(int xi=0;xi<m+n;xi++)
				fin<<btour[xi]<<' ';
			fin<<endl;
			if(eff1>nmax && eff1>nmin)
			{
				//cout<<eff1<<endl;
				//fin<<eff1<<endl;
				nmax=eff1;
				for(xi=0;xi<m+n;xi++)
					btour[xi]=tour[xi];
				//updatepher(tour,eff1);
			}
			if(eff1>evaluate(rt[w-1]))insert(tour,eff1);
			delete tour;
		}
	/*	for(j=0;j<w;j++)
			cout<<evaluate(rt[j])<<endl;
		cout<<endl;
		getch();*/
		for(j=0;j<w;j++)
			updatepher(rt[j],evaluate(rt[j]),j);
		for(j=w-1;j>0;j--)
			rt[j][0]=-1;
		evaporate();
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
	cout<<endl<<" ";
	for(i=0;i<n;i++)
		cout<<" "<<par[i]+1;
	cout<<endl;
	for(i=0;i<m;i++)
		{
		cout<<mac[i]+1<<' ';
		for(int j=0;j<n;j++)
			cout<<pm[mac[i]][par[j]]<<' ';
		cout<<endl;
		}

}

void ants::insert(int *tour,float eff)
{
	for(int i=w-1;i>=0;i--)
	{
		if(eff<evaluate(rt[i]))break;
	}
	int y=i+1;
	for(int j=w-1;j>y;j--)
	{
		if(j>0 && rt[j-1][0]<0)continue;
		for(int k1=0;k1<m+n;k1++)
			rt[j][k1]=rt[j-1][k1];
	}
	for(int k1=0;k1<m+n;k1++)
		rt[y][k1]=tour[k1];
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
/*	for(i=0;i<m;i++)
	{
	for(int j=0;j<n;j++)
		cout<<pm[i][j]<<" ";
	cout<<endl;
	} */
	fin.close();
	km=new long double*[m+n];
	for(i1=0;i1<m+n;i1++)
		km[i1]=new long double[k];
	km1=new long double*[m+n];
	for(i1=0;i1<m+n;i1++)
		km1[i1]=new long double[k];
	for(i1=0;i1<m+n;i1++)
		btour[i1]=0;		
	for(i=0;i<m+n;i++)
		for(j=0;j<k;j++)
			km1[i][j]=km[i][j]=1;
	cout<<"Enter number of ants";
	cin>>an;
	cout<<"Enter the number of best ranked ants out of the above "<<an<<" : ";
	cin>>w;
	rt=new int*[w];
	for(i=0;i<w;i++)
		rt[i]=new int[m+n];
	for(i=0;i<w;i++)
		rt[i][0]=-1;
//	clrscr();
	//reading from files for input complete
}

void main()
{
	ants ant;
//	getch();
	ant.start();
	getch();
}