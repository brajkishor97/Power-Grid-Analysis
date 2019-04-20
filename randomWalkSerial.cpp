#include<cstring>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<cmath>
#include<omp.h>
#include<time.h>
#include<ctime>
#include<stdio.h>

using namespace std;

int main(int argc, char* argv[])
{
  int i,j,total_node;
  int walks,num_threads;
  double *potential;
  if(argc!=2)
  {
    cout<<"Insufficient arguments"<<endl;
    return 1;
  }

  double *rand_array;
  rand_array=(double*)calloc(10000000,sizeof(double));

  for(i=0;i<10000000;i++)
    rand_array[i]=(rand()%1000)/1000.0;

  //cout<<"Enter number of threads: "<<endl;
  //cin>>num_threads;

  cout<<"Enter no of walks: "<<endl;
  cin>>walks;

  //omp_set_num_threads(num_threads);
  double **conductance_array;
  double **probability_array;
  double***neighbor;
  int    **neighbor_id;
  int    *type;
  double *current_value;
  double *node_potential;
  double *updated_potential;
  double *conductance_sum;
  double *reward;
  int    *no_neighbor;

  FILE *myReadFile;
  char typ;
  int node1,node2,node3,k,max=-1;
  j=0;
  double value;
  myReadFile = fopen (argv[1] , "r");

  if(myReadFile!=NULL)
  {
    while(1)
    {
      fscanf(myReadFile, "%c%*d %d %d %lf\n",&typ,&node1,&node2,&value);
      if (typ=='.')
        break;
	    
      node3=(node1>node2)?node1:node2;
      if(node3>max)
      {
        max=node3;
      }
    }
  }
  fclose(myReadFile);

  neighbor=(double***)calloc(max+1,sizeof(double**));
  neighbor_id=(int**)calloc(max+1,sizeof(int*));
  conductance_array=(double**)calloc(max+1,sizeof(double*));
  probability_array=(double**)calloc(max+1,sizeof(double*));
  type=(int*)calloc(max+1,sizeof(int));
  current_value=(double*)calloc(max+1,sizeof(double));
  node_potential=(double*)calloc(max+1,sizeof(double));
  updated_potential=(double*)calloc(max+1,sizeof(double));
  conductance_sum=(double*)calloc(max+1,sizeof(double));
  reward=(double*)calloc(max+1,sizeof(double));
  no_neighbor=(int*)calloc(max+1,sizeof(int));

  for(i=0;i<=max;i++)
  {
    neighbor[i]=(double**)calloc(4,sizeof(double*));
    neighbor_id[i]=(int*)calloc(4,sizeof(int));
    conductance_array[i]=(double*)calloc(4,sizeof(double));
    probability_array[i]=(double*)calloc(4,sizeof(double));
  }

  int no1,no2;

  myReadFile = fopen (argv[1] , "r");
  if(myReadFile!=NULL)
  {
    while(1)
    {
      fscanf(myReadFile, "%c%*d %d %d %lf\n",&typ,&node1,&node2,&value);
      if (typ=='.')
        break;
      no1=no_neighbor[node1];
      no2=no_neighbor[node2];

      if(typ=='R')
      {
        neighbor[node1][no1]=&reward[node2];
        neighbor[node2][no2]=&reward[node1];
        neighbor_id[node1][no1]=node2;
	neighbor_id[node2][no2]=node1;
        conductance_array[node1][no1]=1/value;
        conductance_array[node2][no2]=1/value;
        no_neighbor[node1]++;
        no_neighbor[node2]++;
      }

      else if(typ=='V')
      {
        node_potential[node1]=value;
        type[node1]=1;
      }
      else
      {
        if(node2==0)
          current_value[node1]=value;
        else
          current_value[node2]=-value;
      }
    }
  }

  fclose(myReadFile);

  total_node=max+1;
  cout<<"total node = "<<total_node<<endl;

  for(i=1;i<total_node;i++)
    conductance_sum[i]=conductance_array[i][0]+conductance_array[i][1]+conductance_array[i][2]+conductance_array[i][3];

  for(i=1;i<total_node;i++)
  {
    double temp=0;
    for(j=0;j<no_neighbor[i];j++)
    {
      probability_array[i][j]=temp+conductance_array[i][j]/conductance_sum[i];
      //cerr<<"probability "<<probability_array[i][j]<<endl;
      temp=probability_array[i][j];
    }
  }

  for(i=1;i<total_node;i++)
  {
    if(type[i]==1)
      reward[i]=node_potential[i];
    else
      reward[i]=-current_value[i]/conductance_sum[i];
  }

  int l,iter,tot_iter=0,node_iter;
  long int ind=0;
  double tot_reward,rand_no;
  cerr<<"Started solving: ";
  //double time_st=dsecnd();

  int start_s=clock();

  //#pragma omp parallel for private (i,j,k,l,tot_reward,rand_no,iter,node_iter) firstprivate(walks,total_node) reduction(+:ind)

  for(i=1;i<total_node;i++)
  {
    //node_iter=0;
    for(j=0;j<walks;j++)
    {
      tot_reward=reward[i];
      l=i;
      iter=0;
      while(1)
      {
        if(type[l]==1)
          break;
        rand_no=rand_array[ind%10000000];
        //cerr<<" "<<i<<" neighbor "<<l<<" "<<rand_no<< endl;
        //cerr<<" "<<no_neighbor[l]<<endl;
        for(k=0;k<no_neighbor[l];k++)
        {
          if(rand_no<=probability_array[l][k])
	  {
            break;
          }
        }
        //cerr<<" "<<k<<endl;
        tot_reward=tot_reward+*neighbor[l][k];
        //cout<<"neighbor["<<l<<"]["<<k<<"]="<<*neighbor[l][k]<<endl;
        //cerr<<" "<<iter<<" "<<tot_reward<<endl;
        l=neighbor_id[l][k];
        ++iter;
        ++ind;
      }
      updated_potential[i]=updated_potential[i]+tot_reward;
    }
    node_potential[i]=updated_potential[i]/walks;
    type[i]=1;
    reward[i]=node_potential[i];
  }

  cout<<endl;
  cout<<"Total node traversed: "<<ind<<endl;

  int stop_s=clock();
  cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;
  //cerr<<"Solution"<<endl;

  ofstream arrayData("output.txt",ios::trunc); // File Creation(on C drive)
  for(i=1;i<total_node;i++)
  {
    arrayData<<"V"<<i<<": "<<node_potential[i]<<endl; //Outputs array to txtFile
  }
  return 0;
}
