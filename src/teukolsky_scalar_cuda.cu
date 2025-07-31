#include "teukolsky_scalar_cuda.cuh"
#include "cuda_wrapper.cuh"
#include "pde_cuda_kernel.cuh"



CudaTeukolskyScalarPDE::CudaTeukolskyScalarPDE(Param param_) : param(param_) {
  using namespace Eigen;
  //using namespace Teukolsky;
    
  const Scalar rast_min = param.rast_min;
  const Scalar rast_max = param.rast_max;
  const auto N = param.N;
  const Scalar M = param.M;
  const Scalar a = param.a;
  const auto s = param.s;
  const auto l_max = param.l_max;

  grid_size = N + 1;
  lm_size = (l_max + 1) * (l_max + 1);
    
  const Scalar h = (rast_max - rast_min) / (N - 1);

  // Load coupling mapping info
  psi_lm_map = Teukolsky::make_coupling_info_map(Teukolsky::psi_lm_coupling_info_scalar, l_max);
  dr_psi_lm_map = Teukolsky::make_coupling_info_map(Teukolsky::dr_psi_lm_coupling_info_scalar, l_max);
  drdr_psi_lm_map = Teukolsky::make_coupling_info_map(Teukolsky::drdr_psi_lm_coupling_info_scalar, l_max);
  dt_psi_lm_map = Teukolsky::make_coupling_info_map(Teukolsky::dt_psi_lm_coupling_info_scalar, l_max);

  // Compute the radial coordinate and coupling coefficients
  std::cout << "(CudaTeukolskyScalarPDE) Preparing coupling coefficients" << std::endl;
  
  std::vector<Teukolsky::ComplexVector> coeffs_eigen;
  Teukolsky::prepare_coeffs_scalar(coeffs_eigen, rast_min, rast_max, N, M, a);
  
  coeffs.resize(coeffs_eigen.size());
  for(size_t i = 0; i < coeffs_eigen.size(); ++i){
    coeffs[i].resize(coeffs_eigen[i].size());
    copy_vector(coeffs[i], coeffs_eigen[i]);
  }

  //Q = [N](const Scalar t)->Vector{ return Vector::Zero(N+1); };

  // Prepare buffer for first and second derivatives of psi_lm
  drdr_psi_lm.resize(lm_size * grid_size);
  dr_psi_lm.resize(lm_size * grid_size);


  // std::cout << "(CudaTeukolskyScalarPDE) Initializing CUDA graph" << std::endl;  
  // State x(2 * lm_size * grid_size);   // Placeholder states
  // State dxdt(2 * lm_size * grid_size);
  // cudaGraph_t system_graph = prepare_cuda_graph(x, dxdt);
  // cudaError_t err = cudaGraphInstantiate(&system_graph_exec, system_graph, 0);
  // std::cout << "(cudaGraphInstantiate) err = " << err << std::endl;
  // err = cudaGraphDestroy(system_graph);
  // std::cout << "(cudaGraphDestroy) err = " << err << std::endl;
}


cudaGraph_t CudaTeukolskyScalarPDE::prepare_cuda_graph(const State &x, State &dxdt)
{
  const Scalar rast_min = param.rast_min;
  const Scalar rast_max = param.rast_max;
  const auto N = param.N;
  const Scalar h = (rast_max - rast_min) / (N - 1);
  
  // Prepare CUDA graph for operator()
  cudaGraph_t graph;
  cudaError_t err = cudaGraphCreate(&graph, 0);

  // std::cout << "(GraphCreate) err = " << err << std::endl;

  cudaGraphNode_t copy_time_derivative_node;
  std::vector<cudaGraphNode_t> compute_derivative_nodes(2 * lm_size);
  cudaGraphNode_t empty_barrier_node;
  std::vector<cudaGraphNode_t> compute_dtdt_nodes(lm_size);
  
  // Add computational of derivatives into the graph

  // Copy first order time derivative
  const long long int dt_grid_begin = lm_size * grid_size;
  err = cudaGraphAddMemcpyNode1D(&copy_time_derivative_node, graph, NULL, 0,
				 (void *)thrust::raw_pointer_cast(dxdt.data()),
				 (const void *)(thrust::raw_pointer_cast(x.data()) + dt_grid_begin),
				 dt_grid_begin * sizeof(thrust::complex<double>),
				 cudaMemcpyDeviceToDevice);

  // std::cout << "(cudaGraphAddMemcpyNode1D) err = " << err << std::endl;
  
  // Compute first and second order spatial derivative  
  for(size_t lm = 0; lm < lm_size; ++lm){
    {
      auto arg1 = thrust::raw_pointer_cast(drdr_psi_lm.data() + lm * grid_size);
      auto arg2 = thrust::raw_pointer_cast(x.data() + lm * grid_size);
      int arg3 = grid_size;
      double arg4 = 1 / (h * h);

      void *ptrs[4] = {(void *)&arg1, (void *)&arg2, (void *)&arg3, (void *)&arg4};
      void **ptrs_casted = (void **)ptrs;
    
      const int threadsPerBlock = 512;
      const int numBlocks = (grid_size + threadsPerBlock - 1) / threadsPerBlock;
    
      cudaKernelNodeParams node_params;
      node_params.func = (void *)CUDAKernel::drdr_complex_double_kernel;
      node_params.gridDim = dim3(numBlocks);
      node_params.blockDim = dim3(threadsPerBlock);
      node_params.sharedMemBytes = 0;
      node_params.kernelParams = ptrs_casted;
      node_params.extra = NULL;
    
      err = cudaGraphAddKernelNode(&compute_derivative_nodes[lm], graph, NULL, 0, &node_params);
      // std::cout << "(cudaGraphAddKernelNode) err = " << err << std::endl;
    }
    
    {
      auto arg1 = thrust::raw_pointer_cast(dr_psi_lm.data() + lm * grid_size);
      auto arg2 = thrust::raw_pointer_cast(x.data() + lm * grid_size);
      int arg3 = grid_size;
      double arg4 = 1 / h;

      void *ptrs[4] = {(void *)&arg1, (void *)&arg2, (void *)&arg3, (void *)&arg4};
      void **ptrs_casted = (void **)ptrs;
    
      const int threadsPerBlock = 512;
      const int numBlocks = (grid_size + threadsPerBlock - 1) / threadsPerBlock;
    
      cudaKernelNodeParams node_params;
      node_params.func = (void *)CUDAKernel::dr_complex_double_kernel;
      node_params.gridDim = dim3(numBlocks);
      node_params.blockDim = dim3(threadsPerBlock);
      node_params.sharedMemBytes = 0;
      node_params.kernelParams = ptrs_casted;
      node_params.extra = NULL;
    
      err = cudaGraphAddKernelNode(&compute_derivative_nodes[lm_size + lm], graph, NULL, 0, &node_params);
      // std::cout << "(cudaGraphAddKernelNode) err = " << err << std::endl;
    }
  }

  // Barrier node
  cudaGraphAddEmptyNode(&empty_barrier_node, graph, compute_derivative_nodes.data(), 2 * lm_size);
  // std::cout << "(cudaGraphAddEmptyNode) err = " << err << std::endl;

  // Assign second order time derivatives  
  for(size_t lm = 0; lm < lm_size; ++lm){
    std::vector<void *> args;
      
    for(auto [lm1, idx1] : psi_lm_map[lm]){
      thrust::complex<double> *coeff_ptr = thrust::raw_pointer_cast(coeffs[idx1].data());
      const thrust::complex<double> *var_ptr = thrust::raw_pointer_cast(x.data() + lm1 * grid_size);
      args.push_back(reinterpret_cast<void *>(coeff_ptr));
      args.push_back(const_cast<void *>(reinterpret_cast<const void *>(var_ptr)));
    }
    for(auto [lm1, idx1] : dt_psi_lm_map[lm]){
      thrust::complex<double> *coeff_ptr = thrust::raw_pointer_cast(coeffs[idx1].data());
      const thrust::complex<double> *var_ptr = thrust::raw_pointer_cast(x.data() + (lm_size + lm1) * grid_size);
      args.push_back(reinterpret_cast<void *>(coeff_ptr));
      args.push_back(const_cast<void *>(reinterpret_cast<const void *>(var_ptr)));
    }
    for(auto [lm1, idx1] : dr_psi_lm_map[lm]){
      thrust::complex<double> *coeff_ptr = thrust::raw_pointer_cast(coeffs[idx1].data());
      thrust::complex<double> *var_ptr = thrust::raw_pointer_cast(dr_psi_lm.data() + lm1 * grid_size);
      args.push_back(reinterpret_cast<void *>(coeff_ptr));
      args.push_back(reinterpret_cast<void *>(var_ptr));
    }
    for(auto [lm1, idx1] : drdr_psi_lm_map[lm]){
      thrust::complex<double> *coeff_ptr = thrust::raw_pointer_cast(coeffs[idx1].data());
      thrust::complex<double> *var_ptr = thrust::raw_pointer_cast(drdr_psi_lm.data() + lm1 * grid_size);
      args.push_back(reinterpret_cast<void *>(coeff_ptr));
      args.push_back(reinterpret_cast<void *>(var_ptr));
    }

    auto arg_lhs = thrust::raw_pointer_cast(dxdt.data() + (lm_size + lm) * grid_size);
    int grid_size_store = grid_size;
    // args.push_back(reinterpret_cast<void *>(&arg_lhs));
    // args.push_back(reinterpret_cast<void *>(&grid_size_store));
    
    // The number of terms on the RHS
    // Since the total number of arguments of the kernel is 2n, num_terms is give below
    const size_t num_terms = args.size() / 2;
    
    std::vector<void *> ptrs(args.size() + 2);
    ptrs[0] = (void *)&arg_lhs;
    ptrs[ptrs.size()-1] = (void *)&grid_size_store;
    for(size_t i = 0; i < args.size(); ++i){
      ptrs[i+1] = (void *)&args[i];
    }

    const int threadsPerBlock = 512;
    const int numBlocks = (grid_size + threadsPerBlock - 1) / threadsPerBlock;

    cudaKernelNodeParams node_params;
    node_params.func = (void *)CUDAKernel::assign_lhs_2terms_complex_double_kernels[num_terms];
    node_params.gridDim = dim3(numBlocks);
    node_params.blockDim = dim3(threadsPerBlock);
    node_params.sharedMemBytes = 0;
    node_params.kernelParams = reinterpret_cast<void **>(ptrs.data());
    node_params.extra = NULL;
    
    // dtdt computation should depend on the empty barrier node
    err = cudaGraphAddKernelNode(&compute_dtdt_nodes[lm], graph, &empty_barrier_node, 1, &node_params);
    // std::cout << "(cudaGraphAddKernelNode) err = " << err << std::endl;
  }

  size_t numNodes;
  err = cudaGraphGetNodes(graph, NULL, &numNodes);
  std::cout << "(prepare_cuda_graph) Prepared new CUDA graph with " << numNodes << " nodes" << std::endl;
  // std::cout << "err = " << err << std::endl;

  return graph;
}

void CudaTeukolskyScalarPDE::operator()(const State &x, State &dxdt, const Scalar t)
{
  // cudaGraphExecUpdate(cudaGraphExec_t hGraphExec, cudaGraph_t hGraph, cudaGraphExecUpdateResultInfo* resultInfo )

  // cudaGraphExecUpdateResultInfo update_info;
  // cudaGraph_t new_graph = prepare_cuda_graph(x, dxdt);
  // cudaGraphExecUpdate(system_graph_exec, new_graph, &update_info);
  // cudaGraphLaunch(system_graph_exec, 0);
  // cudaGraphDestroy(new_graph);
  
  const void *x_ptr = reinterpret_cast<const void *>(thrust::raw_pointer_cast(x.data()));
  const void *dxdt_ptr = reinterpret_cast<const void *>(thrust::raw_pointer_cast(dxdt.data()));
  const std::pair<const void *, const void *> key = std::make_pair(x_ptr, dxdt_ptr);
  if(!graph_exec_mapping.contains(key)){
    cudaGraph_t new_graph = prepare_cuda_graph(x, dxdt);
    graph_exec_mapping[key];
    cudaError_t err = cudaGraphInstantiate(&graph_exec_mapping[key], new_graph, 0);
    cudaGraphDestroy(new_graph);
  }
  cudaGraphLaunch(graph_exec_mapping[key], 0);
}
