Recycle BC
===========

UDF that implements a recycle BC that work in parallel.

#include "udf.h"
#include "stdio.h"
#define ID 16              /**** ID of recycle **/
#define nodeid 2         /**** node_ID of recycle **/
#define n_faces 133   /**** face number on recycle **/

real locatx[1000]={0};
real locaty[1000]={0};
real uz[1000]={0};      /***arrays to save data of reycle **/

real x[ND_ND];
real x_in[ND_ND];
int i;
int j;
int k;

DEFINE_PROFILE(v_profile,t_in,m)
{
#if !RP_HOST
	face_t f,f_in;
	cell_t c0;
	Thread *t0,*t;
	Domain *d;

	d=Get_Domain(1);
	t=Lookup_Thread(d,ID);    /****get the face_thread of recycle **/
#endif
/****************************************************************************/  
/****** get data of recycle plane and save in arrays**************/
#if !RP_HOST
	i=0;
    begin_f_loop(f,t) 
	{
		c0=F_C0(f,t);
		t0=THREAD_T0(t);    /**** get the adjacent cell thread ****/
		C_CENTROID(x,c0,t0);
		locatx[i]=x[0];
		locaty[i]=x[1];
		uz[i]=C_W(c0,t0);  /**** W_velocity will be imposed on inlet as condition boundary ****/
		i++;
	}
	end_f_loop(f,t)
#endif
/**************************************************************************/
/******************** sending data to node of inlet ****************/
#if RP_NODE
		if(I_AM_NODE_LAST_P)  /**** node_ID of recycle ****/
		{
			PRF_CSEND_REAL(node_zero,locatx,n_faces,myid);
			PRF_CSEND_REAL(node_zero,locaty,n_faces,myid);
			PRF_CSEND_REAL(node_zero,uz,n_faces,myid);
		}
#endif

/***************************************************************************/
/*************** recieve data from node_recycle*******************/
#if RP_NODE
		if(I_AM_NODE_ZERO_P)  /*** node_ID of inlet  ****/
		{
			PRF_CRECV_REAL(nodeid,locatx,n_faces,nodeid);
			PRF_CRECV_REAL(nodeid,locaty,n_faces,nodeid);
			PRF_CRECV_REAL(nodeid,uz,n_faces,nodeid);
		}
#endif

/***************************************************************************/
/****************map the inlet boundary condtion *****************/
#if !RP_HOST
		j=0;
		begin_f_loop(f_in,t_in)
		{
			F_CENTROID(x_in,f_in,t_in);
			for(k=0;k<n_faces;k++)    /****  look for the target location ****/
				if(fabs(x_in[0]-locatx[k])<=0.00002)    /*** if the X_coordinate of inlet equals the one of inlet **/
					if(fabs(x_in[1]-locaty[k])<=0.00002)    /**** and Y_coordinates equal ***/
					{
						F_PROFILE(f_in,t_in,m)=1000*uz[k]; 
						/***** map inlet velocity with corresponding W_velocity of recycle *****/
						break;   /*** break out the "for loop"  **/
					}
			j++;
		}
		end_f_loop(f_in,t_in)
#endif
}

