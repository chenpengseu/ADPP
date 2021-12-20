%% CRLB
LAtmp = [];
for idx_k = 1:K
	tmp = diag(steerVec(theta(idx_k), posUAV+d_per)...
		.*steerVec(psi, d_per))*s(idx_k)*1j*2*pi*(sind(theta(idx_k))+sind(psi));
	if idx_k == 1
		LAtmp = tmp;
	else
		LAtmp = LAtmp+tmp;
	end
end
LA = A*LAtmp;


for idx_k=1:K
	tmp = 1j*2*pi*s(idx_k)*cosd(theta(idx_k))*exp(1j*2*pi*(d_per*...
		(sind(theta(idx_k))+sind(psi))+posUAV*sind(theta(idx_k)))).*(posUAV+d_per);
	if idx_k == 1
		UPtmp = zeros(length(tmp), K);
	end
	UPtmp(:, idx_k) = tmp;
end
UP = A*UPtmp;

FF = 2/noiseVar*real([LA'*LA, LA'*UP; UP'*LA, UP'*UP]);
rCRLB = sqrt(diag(pinv(FF)));
rCRLB_d(idx_SNR) = rCRLB_d(idx_SNR)+sum(rCRLB(1:N).^2)/N;
rCRLB_theta(idx_SNR) = rCRLB_theta(idx_SNR)+sum(rCRLB(N+1:N+K).^2)/K;
        
%% the proposed method
ang_range = [-60:1:60].';
for idx_test = 1:length(ang_range)
	tmp1 = steerVec(ang_range(idx_test), posUAV);
	tmp2 = steerVec(ang_range(idx_test), posUAV+d_est);
	if idx_test == 1
		AA_T = zeros(length(tmp1), length(ang_range));
		BB_T = zeros(length(tmp2), length(ang_range));
	end

	AA_T(:, idx_test) = tmp1;
	BB_T(:, idx_test) = tmp2;
end
est_t = pinv(kron(AA_T.',eye(N)))*vec(BB_T);
est_T = reshape(est_t, N, N)';
est_T_pinv = inv(est_T);


[V_eig, D_eig] = eig(inv((A*diag(steerVec(psi, d_est)))'*...
	(A*diag(steerVec(psi, d_est)))));
D_vec = real(diag(D_eig));
D_vec(D_vec<0) = 0;
AA = V_eig*diag(D_vec)*V_eig';
	
% with pertubration
cvx_solver sdpt3
cvx_begin sdp quiet
	variable G(N+1, N+1) hermitian;  
	G>=0;
	G(N+1, N+1) == t;
	trace(G) == 1+t;
	for idx = 1:N-1
		sum(diag(G(1:N, 1:N), idx)) == 0;
	end 
	minimize(quad_form((A*diag(steerVec(psi, d_est)))'*r...
		-est_T_pinv*G(1:N, N+1), AA))
cvx_end
x = G(1:N, N+1);
% plot 
sp = abs(x'*steerVec(ang_grid, posUAV));
sp = sp/max(sp);

est_ang = ang_grid(find_peak(sp, K));
err_doa(idx_iter) = sqrt(norm(est_ang-theta)^2/K);

obj_r = [];
FK = [];

if idx_iter == 1
	alpha = .5e-4;
	d_est = zeros(N, 1);
else 
	alpha = alpha*1;
end
alpha2 = .5e-5;

for iter = 1:1e5
	est_s = pinv(A*diag(steerVec(psi, d_est))*(steerVec(est_ang, posUAV).*...
		steerVec(est_ang, d_est))) * r;
	
	sum_v = conj(B)*B.'*diag((steerVec(est_ang, posUAV)...
		.*steerVec(est_ang, d_est).*steerVec(repmat(psi, length(est_ang), 1), d_est))*...
		diag(1j*2*pi*(sind(psi)+sind(est_ang)))*est_s);
	
	v1 = diag(est_s)' * (steerVec(est_ang, d_est).*steerVec(psi * ones(length(est_ang), 1), d_est).*...
		steerVec(est_ang, posUAV))'* sum_v;
	
	v2= diag(1j*2*pi.*(sind(psi)+sind(est_ang)).*est_s)*...
		(steerVec(est_ang, posUAV).*steerVec(psi*ones(length(est_ang), 1), d_est).*steerVec(est_ang, d_est)).'*diag(r'*B.');
	
	d_delta =  sum(2*real(v1-v2)).';
	
	% doa refinement
	d_r = (A*diag(steerVec(psi, d_est))*(steerVec(est_ang, posUAV).*...
		steerVec(est_ang, d_est))*est_s-r)'*A;
	d_theta = 2*real(1j*2*pi*cosd(est_ang).'.*est_s.'.*(d_r*(steerVec(psi*ones(K, 1), d_est).*steerVec(est_ang, d_est+posUAV).*repmat(posUAV+d_est,1,K))));
	d_theta = d_theta.';
	
	est_ang = est_ang-alpha2*d_theta;
	
	d_est = d_est - alpha*d_delta;
	tmp = B.'*diag(steerVec(psi, d_est))*(steerVec(est_ang, posUAV).* steerVec(est_ang, d_est))*est_s;
	obj_r = [obj_r;norm(tmp-r).^2];
	FK = [FK; sqrt(norm(d_per-d_est)^2/N)];
	if iter>100
		tmp = abs((obj_r(end)-obj_r(end-1))/obj_r(end-1));
		if tmp<1e-4
			break;
		end
	end 
end
