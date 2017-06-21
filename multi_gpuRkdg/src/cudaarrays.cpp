#include "../inc/cudaarrays.h"

CCUDAArrays::CCUDAArrays():
area(NULL),		
perimeter(NULL),		
outer_normal_vector(NULL),
mass_coeff(NULL),

vol_bf_value(NULL),	
vol_bdf_value(NULL),
edge_bf_value(NULL),

vol_gauss_weight(NULL),
edge_gauss_weight(NULL),

freedom_rho(NULL),
freedom_rhou(NULL),
freedom_rhov(NULL),
freedom_rhoE(NULL),

convar_rho_vol(NULL),
convar_rhou_vol(NULL),
convar_rhov_vol(NULL),
convar_rhoE_vol(NULL),
convar_rho_edge(NULL),
convar_rhou_edge(NULL),
convar_rhov_edge(NULL),
convar_rhoE_edge(NULL),

rhs_volume_rho(NULL),
rhs_volume_rhou(NULL),
rhs_volume_rhov(NULL),
rhs_volume_rhoE(NULL),
rhs_edge_rho(NULL),
rhs_edge_rhou(NULL),
rhs_edge_rhov(NULL),
rhs_edge_rhoE(NULL),

lfflux_coeff(NULL),

fedge_rho(NULL),	
fedge_rhou(NULL),
fedge_rhov(NULL),
fedge_rhoE(NULL),
gedge_rho(NULL),
gedge_rhou(NULL),
gedge_rhov(NULL),
gedge_rhoE(NULL),

lfflux_rho(NULL),
lfflux_rhou(NULL),
lfflux_rhov(NULL),
lfflux_rhoE(NULL),

freedom_rho_old(NULL),
freedom_rhou_old(NULL),
freedom_rhov_old(NULL),
freedom_rhoE_old(NULL),

ddt(NULL),
residual(NULL),
neighbour(NULL),
sharedEdge(NULL),
triangle_flag(NULL),

_int_pitch(0),
_double_pitch(0)
{

}

void CCUDAArrays::allocateMemory(int num)
{
	size_t width = sizeof(int)*num;

	cudaMallocPitch((void**)&neighbour,     &_int_pitch, width, TRIANGLE_EDGES);
	cudaMallocPitch((void**)&sharedEdge,    &_int_pitch, width, TRIANGLE_EDGES);
	cudaMallocPitch((void**)&triangle_flag, &_int_pitch, width, 1);

	size_t gsize = sizeof(double)*num;
	
	cudaMallocPitch((void**)&area,      &_double_pitch, gsize, 1);
	cudaMallocPitch((void**)&perimeter, &_double_pitch, gsize, 1);

	cudaMallocPitch((void**)&outer_normal_vector, &_double_pitch, gsize, TRIANGLE_EDGES*2); // xy两个分量

	cudaMallocPitch((void**)&mass_coeff, &_double_pitch, gsize, BASIS_FUNCTIONS);
	
	cudaMallocPitch((void**)&vol_bf_value,  &_double_pitch, gsize, VOLUME_GPOINTS*BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&vol_bdf_value, &_double_pitch, gsize, VOLUME_GPOINTS*BASIS_FUNCTIONS*2); // 对x,y的偏导
	cudaMallocPitch((void**)&edge_bf_value, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS*BASIS_FUNCTIONS);
	
	cudaMallocPitch((void**)&vol_gauss_weight,  &_double_pitch, gsize, VOLUME_GPOINTS);
	cudaMallocPitch((void**)&edge_gauss_weight, &_double_pitch, gsize, EDGE_GPOINTS*TRIANGLE_EDGES);


	cudaMallocPitch((void**)&freedom_rho,  &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&freedom_rhou, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&freedom_rhov, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&freedom_rhoE, &_double_pitch, gsize, BASIS_FUNCTIONS);

	cudaMallocPitch((void**)&convar_rho_vol,  &_double_pitch, gsize, VOLUME_GPOINTS);
	cudaMallocPitch((void**)&convar_rhou_vol, &_double_pitch, gsize, VOLUME_GPOINTS);
	cudaMallocPitch((void**)&convar_rhov_vol, &_double_pitch, gsize, VOLUME_GPOINTS);
	cudaMallocPitch((void**)&convar_rhoE_vol, &_double_pitch, gsize, VOLUME_GPOINTS);

	cudaMallocPitch((void**)&convar_rho_edge,  &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&convar_rhou_edge, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&convar_rhov_edge, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&convar_rhoE_edge, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMallocPitch((void**)&rhs_volume_rho,  &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&rhs_volume_rhou, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&rhs_volume_rhov, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&rhs_volume_rhoE, &_double_pitch, gsize, BASIS_FUNCTIONS);
	
	cudaMallocPitch((void**)&rhs_edge_rho,  &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&rhs_edge_rhou, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&rhs_edge_rhov, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&rhs_edge_rhoE, &_double_pitch, gsize, BASIS_FUNCTIONS);

	cudaMallocPitch((void**)&lfflux_coeff, &_double_pitch, gsize, TRIANGLE_EDGES);

	cudaMallocPitch((void**)&fedge_rho,  &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&fedge_rhou, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&fedge_rhov, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&fedge_rhoE, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMallocPitch((void**)&gedge_rho,  &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&gedge_rhou, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&gedge_rhov, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&gedge_rhoE, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMallocPitch((void**)&lfflux_rho,  &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&lfflux_rhou, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&lfflux_rhov, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMallocPitch((void**)&lfflux_rhoE, &_double_pitch, gsize, TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMallocPitch((void**)&freedom_rho_old,  &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&freedom_rhou_old, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&freedom_rhov_old, &_double_pitch, gsize, BASIS_FUNCTIONS);
	cudaMallocPitch((void**)&freedom_rhoE_old, &_double_pitch, gsize, BASIS_FUNCTIONS);

	cudaMalloc((void**)&ddt, sizeof(double));
	cudaMalloc((void**)&residual, sizeof(double)*RESIDUAL_VARS);

	cudaMalloc((void**)&elem_index, width);

	if ( cudaPeekAtLastError()!=cudaSuccess )
		throw CMyException(cudaGetErrorString(cudaPeekAtLastError()));

}

/*
void CCUDAArrays::memsetArrays( int num )
{

	size_t gsize = sizeof(double)*num;

	cudaMemset(freedom_rho,  0, gsize*BASIS_FUNCTIONS);
	cudaMemset(freedom_rhou, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(freedom_rhov, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(freedom_rhoE, 0, gsize*BASIS_FUNCTIONS);

	cudaMemset(convar_rho_vol,  0, gsize*VOLUME_GPOINTS);
	cudaMemset(convar_rhou_vol, 0, gsize*VOLUME_GPOINTS);
	cudaMemset(convar_rhov_vol, 0, gsize*VOLUME_GPOINTS);
	cudaMemset(convar_rhoE_vol, 0, gsize*VOLUME_GPOINTS);

	cudaMemset(convar_rho_edge,  0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(convar_rhou_edge, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(convar_rhov_edge, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(convar_rhoE_edge, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMemset(rhs_volume_rho,  0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rhs_volume_rhou, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rhs_volume_rhov, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rhs_volume_rhoE, 0, gsize*BASIS_FUNCTIONS);

	cudaMemset(rhs_edge_rho,  0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rhs_edge_rhou, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rhs_edge_rhov, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rhs_edge_rhoE, 0, gsize*BASIS_FUNCTIONS);

	cudaMemset(lfflux_coeff, 0, gsize*TRIANGLE_EDGES);

	cudaMemset(fedge_rho,  0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(fedge_rhou, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(fedge_rhov, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(fedge_rhoE, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMemset(gedge_rho,  0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(gedge_rhou, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(gedge_rhov, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(gedge_rhoE, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMemset(lfflux_rho,  0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(lfflux_rhou, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(lfflux_rhov, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);
	cudaMemset(lfflux_rhoE, 0, gsize*TRIANGLE_EDGES*EDGE_GPOINTS);

	cudaMemset(rkdg_rho_t0,  0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rkdg_rhou_t0, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rkdg_rhov_t0, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rkdg_rhoE_t0, 0, gsize*BASIS_FUNCTIONS);

	cudaMemset(rkdg_rho_t1,  0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rkdg_rhou_t1, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rkdg_rhov_t1, 0, gsize*BASIS_FUNCTIONS);
	cudaMemset(rkdg_rhoE_t1, 0, gsize*BASIS_FUNCTIONS);
}
*/

size_t CCUDAArrays::getIntPitch( void ) const
{
	return _int_pitch;
}


size_t CCUDAArrays::getDoublePitch( void ) const
{
	return _double_pitch;
}



CCUDAArrays::~CCUDAArrays()
{
	cudaFree(neighbour);
	cudaFree(sharedEdge);
	cudaFree(triangle_flag);

	cudaFree(area);				
	cudaFree(perimeter);			
	cudaFree(outer_normal_vector);
	cudaFree(mass_coeff);			
	cudaFree(vol_bf_value);		
	cudaFree(vol_bdf_value);		
	cudaFree(edge_bf_value); 		
	cudaFree(vol_gauss_weight);	
	cudaFree(edge_gauss_weight);	

	cudaFree(freedom_rho);	
	cudaFree(freedom_rhou);	
	cudaFree(freedom_rhov);	
	cudaFree(freedom_rhoE);	

	cudaFree(convar_rho_vol);	
	cudaFree(convar_rhou_vol);
	cudaFree(convar_rhov_vol);
	cudaFree(convar_rhoE_vol);
	cudaFree(convar_rho_edge);
	cudaFree(convar_rhou_edge);
	cudaFree(convar_rhov_edge);
	cudaFree(convar_rhoE_edge);

	cudaFree(rhs_volume_rho);	
	cudaFree(rhs_volume_rhou);
	cudaFree(rhs_volume_rhov);
	cudaFree(rhs_volume_rhoE);

	cudaFree(rhs_edge_rho);	
	cudaFree(rhs_edge_rhou);	
	cudaFree(rhs_edge_rhov);	
	cudaFree(rhs_edge_rhoE);	

	cudaFree(lfflux_coeff);	

	cudaFree(fedge_rho);		
	cudaFree(fedge_rhou);		
	cudaFree(fedge_rhov);		
	cudaFree(fedge_rhoE);		

	cudaFree(gedge_rho);		
	cudaFree(gedge_rhou);		
	cudaFree(gedge_rhov);		
	cudaFree(gedge_rhoE);		

	cudaFree(lfflux_rho);		
	cudaFree(lfflux_rhou);	
	cudaFree(lfflux_rhov);	
	cudaFree(lfflux_rhoE);	

	// rkdg 中间量
	cudaFree(freedom_rho_old);	
	cudaFree(freedom_rhou_old);	
	cudaFree(freedom_rhov_old);	
	cudaFree(freedom_rhoE_old);	

	cudaFree(ddt);
	cudaFree(residual);
}