function Xn = ecef2enu(X, P)
R = Get_Rotation_Matrix(P);
Xn = R'*(X-P*ones(1,size(X,2)));
end