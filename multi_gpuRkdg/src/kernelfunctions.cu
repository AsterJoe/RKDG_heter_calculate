#include "../inc/kernelfunctions.cuh"

__global__ void kernel_getTimeStep(
	int tnum, double gamma, double cfl, double *ddt,

	const double *freedom_rho,   const double *freedom_rhou, 
	const double *freedom_rhov,  const double *freedom_rhoE,

	const double *perimeter, const double *area
	)
{
	extern __shared__ double mx[];

	int tid = threadIdx.x;
	double tmx(0);

	for ( int index=tid; index<tnum; index+=blockDim.x )
	{
		double rho = freedom_rho[index];
		double u   = freedom_rhou[index] / rho;
		double v   = freedom_rhov[index] / rho;
		double p   = (gamma-1)*(freedom_rhoE[index]-0.5*rho*(u*u+v*v));
		double a   = sqrt(gamma*p/rho);

		double dx  = (sqrt(u*u+v*v)+a)*perimeter[index]/area[index];

		if ( dx>tmx )
		{
			tmx = dx;
		}
	}

	mx[tid] = tmx;
	__syncthreads();

	// 选出最大的值
	for ( int index=blockDim.x/2; index>64; index/=2 )
	{
		if ( tid<index )
		{
			if ( mx[tid+index]>mx[tid] )
			{
				mx[tid] = mx[tid+index];
			}
		}
		__syncthreads();
	}

	// 剩下的32个线程
	if ( tid<32 )
	{
		volatile double *smx = mx;
		if ( smx[tid+32]>smx[tid] )
		{
			smx[tid] = smx[tid+32];
		}
		if ( smx[tid+16]>smx[tid] )
		{
			smx[tid] = smx[tid+16];
		}
		if ( smx[tid+8]>smx[tid] )
		{
			smx[tid] = smx[tid+8];
		}
		if ( smx[tid+4]>smx[tid] )
		{
			smx[tid] = smx[tid+4];
		}
		if ( smx[tid+2]>smx[tid] )
		{
			smx[tid] = smx[tid+2];
		}
		if ( smx[tid+1]>smx[tid] )
		{
			smx[tid] = smx[tid+1];
		}

		if ( tid==0 )
			ddt[0] = cfl / smx[0];
	}
}


__global__ void kernel_calculateConVars(int tnum, int double_pitch,

	const double *freedom_rho,  const double *freedom_rhou,
	const double *freedom_rhov, const double *freedom_rhoE,

	double *convar_rho_vol,  double *convar_rhou_vol,
	double *convar_rhov_vol, double *convar_rhoE_vol,

	double *convar_rho_edge,   double *convar_rhou_edge, 
	double *convar_rhov_edge,  double *convar_rhoE_edge, 

	const double *vol_bf_value, const double *edge_bf_value
	)
{
	int tid = threadIdx.x + blockDim.x*blockIdx.x;
	double vc[4], ec[4], bfv;
	extern __shared__ double freedom[];

	if ( tid<tnum )
	{
		// 单元内的守恒量
		for ( int i=0, base_index=tid; i<TRIANGLE_EDGES*EDGE_GPOINTS; ++i, base_index+=double_pitch )
		{
			if ( i<VOLUME_GPOINTS )
			{
				vc[0] = 0; vc[1] = 0; vc[2] = 0; vc[3] = 0;
			}

			ec[0] = 0; ec[1] = 0; ec[2] = 0; ec[3] = 0;


			for ( int j=0, index=tid; j<BASIS_FUNCTIONS; ++j, index+=double_pitch )
			{

				freedom[threadIdx.x]              = freedom_rho[index];
				freedom[threadIdx.x+blockDim.x]   = freedom_rhou[index];
				freedom[threadIdx.x+2*blockDim.x] = freedom_rhov[index];
				freedom[threadIdx.x+3*blockDim.x] = freedom_rhoE[index];

				if ( i<VOLUME_GPOINTS )
				{
					bfv = vol_bf_value[base_index+j*VOLUME_GPOINTS*double_pitch];
					
					vc[0] += bfv*freedom[threadIdx.x]             ;
					vc[1] += bfv*freedom[threadIdx.x+blockDim.x]  ;
					vc[2] += bfv*freedom[threadIdx.x+2*blockDim.x];
					vc[3] += bfv*freedom[threadIdx.x+3*blockDim.x];
				}

				bfv = edge_bf_value[base_index+j*TRIANGLE_EDGES*EDGE_GPOINTS*double_pitch];

				ec[0] += bfv*freedom[threadIdx.x]             ;
				ec[1] += bfv*freedom[threadIdx.x+blockDim.x]  ;
				ec[2] += bfv*freedom[threadIdx.x+2*blockDim.x];
				ec[3] += bfv*freedom[threadIdx.x+3*blockDim.x];

			}

			if ( i<VOLUME_GPOINTS )
			{
				convar_rho_vol[base_index]  = vc[0];
				convar_rhou_vol[base_index] = vc[1];
				convar_rhov_vol[base_index] = vc[2];
				convar_rhoE_vol[base_index] = vc[3];
			}

			convar_rho_edge[base_index]  = ec[0];
			convar_rhou_edge[base_index] = ec[1];
			convar_rhov_edge[base_index] = ec[2];
			convar_rhoE_edge[base_index] = ec[3];
		}
	}
}

__global__ void kernel_boundaryCondition(
	int tnum, int num, int double_pitch, 

	double rho, double rhou, double rhov, double rhoE,

	double *convar_rho_edge, double *convar_rhou_edge, 
	double *convar_rhov_edge, double *convar_rhoE_edge,

	double *freedom_rho, double *freedom_rhou, 
	double *freedom_rhov, double *freedom_rhoE,

	const int *neighbour, const int *sharedEdge, 

	const int *triangle_flag, const double *outer_normal
	)
{
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	int ne, se, rindex;
	int i, j, index;
	double nx, ny;

	if ( tid<num-tnum )
	{
		index = tid + tnum;
		switch ( triangle_flag[index] )
		{
		case CELL_FARFIELD:

			freedom_rho[index]  = rho;
			freedom_rhou[index] = rhou;
			freedom_rhov[index] = rhov;
			freedom_rhoE[index] = rhoE;


			for ( i=0; i<TRIANGLE_EDGES*EDGE_GPOINTS; ++i, index+=double_pitch )
			{
				convar_rho_edge[index]  = rho;
				convar_rhou_edge[index] = rhou;
				convar_rhov_edge[index] = rhov;
				convar_rhoE_edge[index] = rhoE;
				
			}
			
			break;

		case CELL_OUTFLOW:

			ne = neighbour[index];
			se = sharedEdge[index];

			freedom_rho[index]  = freedom_rho[ne] ;
			freedom_rhou[index] = freedom_rhou[ne];
			freedom_rhov[index] = freedom_rhov[ne];
			freedom_rhoE[index] = freedom_rhoE[ne];


			for ( i=0; i<TRIANGLE_EDGES; ++i )
			{
				index = tid+tnum + i*EDGE_GPOINTS*double_pitch;
				rindex = ne+(EDGE_GPOINTS-1)*double_pitch+se*EDGE_GPOINTS*double_pitch;
					
				for ( j=0; j<EDGE_GPOINTS; ++j,index+=double_pitch, rindex-=double_pitch )
				{

					convar_rho_edge[index]  = convar_rho_edge[rindex];
					convar_rhou_edge[index] = convar_rhou_edge[rindex];
					convar_rhov_edge[index] = convar_rhov_edge[rindex];
					convar_rhoE_edge[index] = convar_rhoE_edge[rindex];
				}
				
			}
			
			break;

		case CELL_REFLECTION:

			ne = neighbour[index];
			se = sharedEdge[index];

			nx = outer_normal[ne+2*se*double_pitch];
			ny = outer_normal[ne+(2*se+1)*double_pitch];

			freedom_rho[index]  = freedom_rho[ne] ;
			freedom_rhoE[index] = freedom_rhoE[ne];

			freedom_rhou[index] = freedom_rhou[ne]*(ny*ny-nx*nx)
								- 2.0*nx*ny*freedom_rhov[ne];
			freedom_rhov[index] = freedom_rhov[ne]*(nx*nx-ny*ny)
								- 2.0*nx*ny*freedom_rhou[ne];


			for ( i=0; i<TRIANGLE_EDGES; ++i )
			{
				index = tid+tnum+i*EDGE_GPOINTS*double_pitch;
				rindex = ne+(EDGE_GPOINTS-1)*double_pitch+se*EDGE_GPOINTS*double_pitch;
				
				for ( j=0; j<EDGE_GPOINTS; ++j, index+=double_pitch, rindex-=double_pitch )
				{

					convar_rho_edge[index]  = convar_rho_edge[rindex];
					convar_rhoE_edge[index] = convar_rhoE_edge[rindex];

					// 速度反射
					convar_rhou_edge[index] = convar_rhou_edge[rindex]*(ny*ny-nx*nx)
											- 2.0*nx*ny*convar_rhov_edge[rindex];

					convar_rhov_edge[index] = convar_rhov_edge[rindex]*(nx*nx-ny*ny)
											- 2.0*nx*ny*convar_rhou_edge[rindex];
				}
			}
			break;

		default:
			break;
		}
	}
	
}


__global__ void kernel_calculateVolumeRHS(
	int tnum, int double_pitch, double gamma, 

	const double *convar_rho_vol,const double *convar_rhou_vol, 
	const double *convar_rhov_vol, const double *convar_rhoE_vol,

	double *rhs_volume_rho, double *rhs_volume_rhou, 
	double *rhs_volume_rhov, double *rhs_volume_rhoE,

	const double *vol_gauss_weight, const double *vol_bdf_value
	)
{
	int tid = threadIdx.x + blockDim.x*blockIdx.x;

	int index, i, j;
	double rhs[4], f[4];

	extern __shared__ double sw[];

	if ( tid<tnum )
	{	
		for ( i=0, index=tid; i<VOLUME_GPOINTS; ++i, index+=double_pitch )
			sw[threadIdx.x+i*blockDim.x] = vol_gauss_weight[index]; 

		for ( i=0; i<BASIS_FUNCTIONS; ++i )
		{
			
			rhs[0] = 0; rhs[1] = 0; rhs[2] = 0; rhs[3] = 0;

			for ( j=0, index=tid; j<VOLUME_GPOINTS; ++j, index+=double_pitch )
			{

				double rho  = convar_rho_vol[index];
				double u    = convar_rhou_vol[index]/rho;
				double v    = convar_rhov_vol[index]/rho;
				double rhoE = convar_rhoE_vol[index];

				double w   = sw[threadIdx.x+j*blockDim.x];
				double dvx = vol_bdf_value[index+2*i*VOLUME_GPOINTS*double_pitch];
				double dvy = vol_bdf_value[index+(2*i+1)*VOLUME_GPOINTS*double_pitch];

				double p = (gamma-1)*(rhoE-0.5*rho*(u*u+v*v));

				f[0] = rho*u;
				f[1] = rho*u*u + p;
				f[2] = rho*u*v;
				f[3] = u*(rhoE+p);

				rhs[0] += w*f[0]*dvx;
				rhs[1] += w*f[1]*dvx;
				rhs[2] += w*f[2]*dvx;
				rhs[3] += w*f[3]*dvx;


				f[0] = rho*v;
				f[1] = f[2];
				f[2] = rho*v*v + p;
				f[3] = v*(rhoE+p);
				

				rhs[0] += w*f[0]*dvy;
				rhs[1] += w*f[1]*dvy;
				rhs[2] += w*f[2]*dvy;
				rhs[3] += w*f[3]*dvy;
			}
			index = tid + i*double_pitch;

			rhs_volume_rho[index]  = rhs[0];
			rhs_volume_rhou[index] = rhs[1];
			rhs_volume_rhov[index] = rhs[2];
			rhs_volume_rhoE[index] = rhs[3];
		}
	}
}


__global__ void kernel_calculateLFCoeff(
	int tnum, int int_pitch, int double_pitch,  double gamma, 
	
	const double *outer_normal_vector,  const int *neighbour, 

	const double* freedom_rho, const double *freedom_rhou, 
	const double *freedom_rhov, const double *freedom_rhoE, 

	double *lfflux_coeff
	)
{
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	if ( tid<tnum )
	{

		for ( int i=0, index=tid; i<TRIANGLE_EDGES; ++i, index+=double_pitch )
		{

			int ne = neighbour[tid+i*int_pitch];
			
			double nx = outer_normal_vector[tid+2*i*double_pitch];
			double ny = outer_normal_vector[tid+(2*i+1)*double_pitch];

			double rho  = 0.5*(freedom_rho[tid] +freedom_rho[ne]);
			double u    = 0.5*(freedom_rhou[tid]+freedom_rhou[ne])/rho;
			double v    = 0.5*(freedom_rhov[tid]+freedom_rhov[ne])/rho;
			double rhoE = 0.5*(freedom_rhoE[tid]+freedom_rhoE[ne]);

			double p  = (gamma-1)*(rhoE-0.5*rho*(u*u+v*v));
			double ma = sqrt(gamma*p/rho);

			lfflux_coeff[index] = fabs(u*nx+v*ny) + ma;
		}
	}
	
}

__global__ void kernel_calculateEdgeFG(
	int tnum, int num, int double_pitch, double gamma, 

	const double *convar_rho_edge, const double *convar_rhou_edge,
	const double *convar_rhov_edge, const double *convar_rhoE_edge,

	double *fedge_rho, double *fedge_rhou,
	double *fedge_rhov, double *fedge_rhoE,

	double *gedge_rho, double *gedge_rhou,
	double *gedge_rhov, double *gedge_rhoE
	)
{
	int tid = threadIdx.x + blockDim.x*blockIdx.x;

	if ( tid<num )
	{

		for ( int i=0, index=tid; i<TRIANGLE_EDGES*EDGE_GPOINTS; ++i,index+=double_pitch )
		{

			double rho  = convar_rho_edge[index];
			double u    = convar_rhou_edge[index] / rho;
			double v    = convar_rhov_edge[index] / rho;
			double rhoE = convar_rhoE_edge[index];

			double p = (gamma-1)*(rhoE-0.5*rho*(u*u+v*v));

			fedge_rho[index]  = rho*u;
			fedge_rhou[index] = rho*u*u + p;
			fedge_rhov[index] = rho*u*v;
			fedge_rhoE[index] = u*(rhoE+p);

			gedge_rho[index]  = rho*v;
			gedge_rhou[index] = rho*v*u;
			gedge_rhov[index] = rho*v*v + p;
			gedge_rhoE[index] = v*(rhoE+p);

		}
		
	}
}

/*
__global__ void kernel_calculateFlux(
	int tnum, int int_pitch, int double_pitch, 
	const int *neighbour, const int *sharedEdge,

	const double *convar_rho_edge, const double *convar_rhou_edge,
	const double *convar_rhov_edge, const double *convar_rhoE_edge,

	const double *fedge_rho, const double *fedge_rhou,
	const double *fedge_rhov, const double *fedge_rhoE,

	const double *gedge_rho, const double *gedge_rhou,
	const double *gedge_rhov, const double *gedge_rhoE,

	const double *outer_normal_vector, const double *lfflux_coeff,

	double *lfflux_rho, double *lfflux_rhou,
	double *lfflux_rhov, double *lfflux_rhoE
	)
{
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	if ( tid<tnum )
	{
		for ( int i=0; i<TRIANGLE_EDGES; ++i )
		{
			int rindex = tid + i*int_pitch;
			
			int ne = neighbour[rindex];
			int se = sharedEdge[rindex];
			
			double lfcoeff = lfflux_coeff[tid + i*double_pitch];
			
			double nx = outer_normal_vector[tid+2*i*double_pitch];
			double ny = outer_normal_vector[tid+(2*i+1)*double_pitch];

			int lindex = tid + i*EDGE_GPOINTS*double_pitch;
			
			rindex = ne + (EDGE_GPOINTS-1)*double_pitch + se*EDGE_GPOINTS*double_pitch;
			
			for ( int j=0; j<EDGE_GPOINTS; ++j, lindex+=double_pitch, rindex-=double_pitch )
			{

				lfflux_rho[lindex] = 0.5*(
					  (fedge_rho[lindex]+fedge_rho[rindex])*nx
					+ (gedge_rho[lindex]+gedge_rho[rindex])*ny
					- lfcoeff*(convar_rho_edge[rindex]-convar_rho_edge[lindex])
					);

				lfflux_rhou[lindex] = 0.5*(
					(fedge_rhou[lindex]+fedge_rhou[rindex])*nx
					+ (gedge_rhou[lindex]+gedge_rhou[rindex])*ny
					- lfcoeff*(convar_rhou_edge[rindex]-convar_rhou_edge[lindex])
					);

				lfflux_rhov[lindex] = 0.5*(
					(fedge_rhov[lindex]+fedge_rhov[rindex])*nx
					+ (gedge_rhov[lindex]+gedge_rhov[rindex])*ny
					- lfcoeff*(convar_rhov_edge[rindex]-convar_rhov_edge[lindex])
					);

				lfflux_rhoE[lindex] = 0.5*(
					(fedge_rhoE[lindex]+fedge_rhoE[rindex])*nx
					+ (gedge_rhoE[lindex]+gedge_rhoE[rindex])*ny
					- lfcoeff*(convar_rhoE_edge[rindex]-convar_rhoE_edge[lindex])
					);
			}
			
		}
		
	}
	
}
*/

__global__ void kernel_calculateFlux(
	int tnum, int int_pitch, int double_pitch, 
	const int *neighbour, const int *sharedEdge,

	const double *convar_rho_edge, const double *convar_rhou_edge,
//	const double *convar_rhov_edge, const double *convar_rhoE_edge,

	const double *fedge_rho, const double *fedge_rhou,
//	const double *fedge_rhov, const double *fedge_rhoE,

	const double *gedge_rho, const double *gedge_rhou,
//	const double *gedge_rhov, const double *gedge_rhoE,

	const double *outer_normal_vector, const double *lfflux_coeff,

	double *lfflux_rho, double *lfflux_rhou
//	double *lfflux_rhov, double *lfflux_rhoE
	)
{
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	if ( tid<tnum )
	{
		for ( int i=0; i<TRIANGLE_EDGES; ++i )
		{
			int rindex = tid + i*int_pitch;
			
			int ne = neighbour[rindex];
			int se = sharedEdge[rindex];
			
			double lfcoeff = lfflux_coeff[tid + i*double_pitch];
			
			double nx = outer_normal_vector[tid+2*i*double_pitch];
			double ny = outer_normal_vector[tid+(2*i+1)*double_pitch];

			int lindex = tid + i*EDGE_GPOINTS*double_pitch;
			
			rindex = ne + (EDGE_GPOINTS-1)*double_pitch + se*EDGE_GPOINTS*double_pitch;
			
			for ( int j=0; j<EDGE_GPOINTS; ++j, lindex+=double_pitch, rindex-=double_pitch )
			{

				lfflux_rho[lindex] = 0.5*(
					  (fedge_rho[lindex]+fedge_rho[rindex])*nx
					+ (gedge_rho[lindex]+gedge_rho[rindex])*ny
					- lfcoeff*(convar_rho_edge[rindex]-convar_rho_edge[lindex])
					);

				lfflux_rhou[lindex] = 0.5*(
					(fedge_rhou[lindex]+fedge_rhou[rindex])*nx
					+ (gedge_rhou[lindex]+gedge_rhou[rindex])*ny
					- lfcoeff*(convar_rhou_edge[rindex]-convar_rhou_edge[lindex])
					);
/*
				lfflux_rhov[lindex] = 0.5*(
					(fedge_rhov[lindex]+fedge_rhov[rindex])*nx
					+ (gedge_rhov[lindex]+gedge_rhov[rindex])*ny
					- lfcoeff*(convar_rhov_edge[rindex]-convar_rhov_edge[lindex])
					);

				lfflux_rhoE[lindex] = 0.5*(
					(fedge_rhoE[lindex]+fedge_rhoE[rindex])*nx
					+ (gedge_rhoE[lindex]+gedge_rhoE[rindex])*ny
					- lfcoeff*(convar_rhoE_edge[rindex]-convar_rhoE_edge[lindex])
					);
*/
			}
			
		}
		
	}
	
}

__global__ void kernel_calculateEdgeRHS(
	int tnum, int double_pitch, 
	const double *edge_gauss_weight, const double *edge_bf_value,

	const double *lfflux_rho,  const double *lfflux_rhou, 
	const double *lfflux_rhov, const double *lfflux_rhoE, 

	double *rhs_edge_rho,  double *rhs_edge_rhou, 
	double *rhs_edge_rhov, double *rhs_edge_rhoE,

	const double *rhs_volume_rho,  const double *rhs_volume_rhou, 
	const double *rhs_volume_rhov, const double *rhs_volume_rhoE
	)
{
	int tid = threadIdx.x + blockDim.x*blockIdx.x;

	double rhs[4],w, bfv;
	int i, j, index;

	extern __shared__ double sw[];

	if ( tid<tnum )
	{
		// 将高斯积分权重读取出来
		for ( i=0, index=tid; i<TRIANGLE_EDGES*EDGE_GPOINTS; ++i,index+=double_pitch )
			sw[threadIdx.x+i*blockDim.x] = edge_gauss_weight[index];

		for ( i=0; i<BASIS_FUNCTIONS; ++i )
		{
			rhs[0] = 0; rhs[1] = 0; rhs[2] = 0; rhs[3] = 0;

			for ( j=0,index=tid; j<TRIANGLE_EDGES*EDGE_GPOINTS; ++j,index+=double_pitch )
			{

				w     = sw[threadIdx.x+j*blockDim.x];
				bfv   = edge_bf_value[index+i*TRIANGLE_EDGES*EDGE_GPOINTS*double_pitch];

				rhs[0] += w*lfflux_rho[index]*bfv;
				rhs[1] += w*lfflux_rhou[index]*bfv;
				rhs[2] += w*lfflux_rhov[index]*bfv;
				rhs[3] += w*lfflux_rhoE[index]*bfv;
			}
			
			index = tid + i*double_pitch;

			rhs_edge_rho[index]  = rhs_volume_rho[index]  - rhs[0];
			rhs_edge_rhou[index] = rhs_volume_rhou[index] - rhs[1];
			rhs_edge_rhov[index] = rhs_volume_rhov[index] - rhs[2];
			rhs_edge_rhoE[index] = rhs_volume_rhoE[index] - rhs[3];
		}
		
	}
	
}


__global__ void kernel_rkdgStepOne(
	int tnum, int double_pitch, double dt, const double *mass_coeff,

	double *freedom_rho,  double *freedom_rhou, 
	double *freedom_rhov, double *freedom_rhoE, 

	const double *rhs_edge_rho,  const double *rhs_edge_rhou,
	const double *rhs_edge_rhov, const double *rhs_edge_rhoE
	)
{
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	if ( tid<tnum )
	{

		for ( int i=0, index=tid; i<BASIS_FUNCTIONS; ++i,index+=double_pitch )
		{

			double mass  = mass_coeff[index];

			freedom_rho[index]  += dt*(rhs_edge_rho[index])/mass;
			freedom_rhou[index] += dt*(rhs_edge_rhou[index])/mass;
			freedom_rhov[index] += dt*(rhs_edge_rhov[index])/mass;
			freedom_rhoE[index] += dt*(rhs_edge_rhoE[index])/mass;
		}
		
	}
	
}


__global__ void kernel_rkdgStepTwo(
	int tnum, int double_pitch, double dt, const double *mass_coeff,
	
	double *freedom_rho,  double *freedom_rhou, 
	double *freedom_rhov, double *freedom_rhoE, 

	const double *rhs_edge_rho,  const double *rhs_edge_rhou,
	const double *rhs_edge_rhov, const double *rhs_edge_rhoE,

	const double *freedom_rho_old,  const double *freedom_rhou_old, 
	const double *freedom_rhov_old, const double *freedom_rhoE_old
	)
{
	int tid = threadIdx.x + blockDim.x*blockIdx.x;

	if ( tid<tnum )
	{

		for ( int i=0,index=tid; i<BASIS_FUNCTIONS; ++i,index+=double_pitch )
		{

			double mass = mass_coeff[index];

			freedom_rho[index]  = 0.75*freedom_rho_old[index] + 0.25*freedom_rho[index]
								+ 0.25*dt*(rhs_edge_rho[index] )/mass;

			freedom_rhou[index] = 0.75*freedom_rhou_old[index] + 0.25*freedom_rhou[index]
								+ 0.25*dt*(rhs_edge_rhou[index] )/mass;

			freedom_rhov[index] = 0.75*freedom_rhov_old[index] + 0.25*freedom_rhov[index]
								+ 0.25*dt*(rhs_edge_rhov[index] )/mass;

			freedom_rhoE[index] = 0.75*freedom_rhoE_old[index] + 0.25*freedom_rhoE[index]
								+ 0.25*dt*(rhs_edge_rhoE[index] )/mass;
		}
		
	}
	
}


__global__ void kernel_rkdgStepThree(
	int tnum, int double_pitch, double dt, const double *mass_coeff,

	double *freedom_rho,  double *freedom_rhou, 
	double *freedom_rhov, double *freedom_rhoE, 

	const double *rhs_edge_rho,  const double *rhs_edge_rhou,
	const double *rhs_edge_rhov, const double *rhs_edge_rhoE,

	const double *freedom_rho_old,  const double *freedom_rhou_old, 
	const double *freedom_rhov_old, const double *freedom_rhoE_old
	)
{
	int tid = threadIdx.x + blockDim.x*blockIdx.x;

	if ( tid<tnum )
	{

		for ( int i=0,index=tid; i<BASIS_FUNCTIONS; ++i,index+=double_pitch )
		{
			double mass = mass_coeff[index];

			freedom_rho[index]  = 1.0/3.0*freedom_rho_old[index] + 2.0/3.0*freedom_rho[index]
						+ 2.0/3.0*dt*(rhs_edge_rho[index] )/mass;

			freedom_rhou[index] = 1.0/3.0*freedom_rhou_old[index] + 2.0/3.0*freedom_rhou[index]
						+ 2.0/3.0*dt*(rhs_edge_rhou[index] )/mass;

			freedom_rhov[index] = 1.0/3.0*freedom_rhov_old[index] + 2.0/3.0*freedom_rhov[index]
						+ 2.0/3.0*dt*(rhs_edge_rhov[index] )/mass;

			freedom_rhoE[index] = 1.0/3.0*freedom_rhoE_old[index] + 2.0/3.0*freedom_rhoE[index]
						+ 2.0/3.0*dt*(rhs_edge_rhoE[index] )/mass;
		}

	}
}


__global__ void kernel_calculateResidual(
	int tnum, 
	const double *freedom_rho,     const double *freedom_rhoE,
	const double *freedom_rho_old, const double *freedom_rhoE_old,
	double *residual
	)
{
	extern __shared__ double res[];

	int tid = threadIdx.x;

	double rhomax(0), Emax(0), ds;
	int index;


	for ( index=tid; index<tnum; index+=blockDim.x )
	{
		ds = fabs(freedom_rho[index]-freedom_rho_old[index]);
		if ( rhomax<ds )
		{
			rhomax = ds;
		}
		
		ds = fabs(freedom_rhoE[index]-freedom_rhoE_old[index]);
		if ( Emax<ds )
		{
			Emax = ds;
		}
		
	}

	res[2*tid]   = rhomax;
	res[2*tid+1] = Emax;
	__syncthreads();

	// 选出最大的值
	for ( index=blockDim.x/2; index>64; index/=2 )
	{
		if ( tid<index )
		{
			if ( res[2*(tid+index)]>res[2*tid] )
			{
				res[2*tid] = res[2*(tid+index)];
			}
			if ( res[2*(tid+index)+1]>res[2*tid+1] )
			{
				res[2*tid+1] = res[2*(tid+index)+1];
			}
		}
		__syncthreads();
	}

	// 剩下的32个线程
	if ( tid<32 )
	{
		volatile double *sres = res;
		if ( sres[2*(tid+32)]>sres[2*tid] )
		{
			sres[2*tid] = sres[2*(tid+32)];
		}
		if ( sres[2*(tid+32)+1]>sres[2*tid+1] )
		{
			sres[2*tid+1] = sres[2*(tid+32)+1];
		}

		if ( sres[2*(tid+16)]>sres[2*tid] )
		{
			sres[2*tid] = sres[2*(tid+16)];
		}
		if ( sres[2*(tid+16)+1]>sres[2*tid+1] )
		{
			sres[2*tid+1] = sres[2*(tid+16)+1];
		}

		if ( sres[2*(tid+8)]>sres[2*tid] )
		{
			sres[2*tid] = sres[2*(tid+8)];
		}
		if ( sres[2*(tid+8)+1]>sres[2*tid+1] )
		{
			sres[2*tid+1] = sres[2*(tid+8)+1];
		}

		if ( sres[2*(tid+4)]>sres[2*tid] )
		{
			sres[2*tid] = sres[2*(tid+4)];
		}
		if ( sres[2*(tid+4)+1]>sres[2*tid+1] )
		{
			sres[2*tid+1] = sres[2*(tid+4)+1];
		}

		if ( sres[2*(tid+2)]>sres[2*tid] )
		{
			sres[2*tid] = sres[2*(tid+2)];
		}
		if ( sres[2*(tid+2)+1]>sres[2*tid+1] )
		{
			sres[2*tid+1] = sres[2*(tid+2)+1];
		}


		if ( sres[2*(tid+1)]>sres[2*tid] )
		{
			sres[2*tid] = sres[2*(tid+1)];
		}
		if ( sres[2*(tid+1)+1]>sres[2*tid+1] )
		{
			sres[2*tid+1] = sres[2*(tid+1)+1];
		}
		if ( tid==0 )
		{
			residual[0] = sres[0];
			residual[1] = sres[1];
		}
	}
}