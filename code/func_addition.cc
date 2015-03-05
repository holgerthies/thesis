alg_func_ptr<PARAM, RESULT> sum(new const alg_func<PARAM, RESULT>(
	[alg1, alg2](const PARAM& x) -> RESULT 
		{return (*alg1)(x) + (*alg2)(x);}
));