function L = scanPolyCurve(m, S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Call: L = scanPolyCurve(m, S)
% Description :  scanPolyCurve creates an image L of size S with values
%                0 on the boundary of the curve m
%                1 outside of the curve m
%                2 inside the curve m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m(:,1) = max(0, min(S(1), m(:,1))) ;
m(:,2) = max(0, min(S(2), m(:,2))) ;


I = ones(S) ;
for k=1:size(m,1)-1,
    d = norm(m(k+1, :)-m(k,:)) ;
    i1 = round(m(k,1) + [0:2*d] *(m(k+1, 1)-m(k,1))/(2*d)) ;
    i2 = round(m(k,2) + [0:2*d] * (m(k+1, 2)-m(k,2))/(2*d)) ;
   
    for l=1:numel(i1),
         if (i1(l)<=0)
	    i1(l) = 1 ;
         elseif (i1(l)>size(I,1))
	    i1(l) = size(I,1)
         elseif (i2(l)<=0)
            i2(l) = 1;
         elseif (i2(l)>size(I,2))
            i2(l) = size(I,2);
         end
         I(i1(l), i2(l)) = 0 ;
    end
end
k = size(m,1);
d = norm(m(1, :)-m(k,:)) + 1e-10 ;
i1 = round(m(k,1) + (0:2*d) * (m(1, 1)-m(k,1))/(2*d)) ;
i2 = round(m(k,2) + (0:2*d) * (m(1, 2)-m(k,2))/(2*d)) ;
for l=1:numel(i1),
         if (i1(l)<=0)
	    i1(l) = 1 ;
         elseif (i1(l)>size(I,1))
	    i1(l) = size(I,1)
         elseif (i2(l)<=0)
            i2(l) = 1;
         elseif (i2(l)>size(I,2))
            i2(l) = size(I,2);
         end
         I(i1(l), i2(l)) = 0 ;
end



L = zeros(S) ;

L(1, :) = 1 ;
L(:,1) = 1 ;

prev = 1 ;
for r = 2:S(1),
    L(r, 1) = 1 ;
    for c = 2:S(2),
        if (I(r,c) ~= 0)
            if (L(r, c-1) ~= 0)
                L(r,c) = L(r,c-1) ;
            else if (L(r-1, c) ~= 0)
                    L(r,c) = L(r-1, c) ;
                    %prev = L(r,c) ;
                else
                    L(r,c) = prev + 1 ;
                    prev = prev+ 1 ;
                end;
            end ;
        end ;
    end ;
end ;

N = max(max(L)) ;
J = 1:N ;
for r = S(1)-1:-1:1,
    for c = S(2)-1:-1:1,
        if (I(r,c) ~= 0)
            if (L(r, c+1) ~= 0 & J(L(r,c)) ~= J(L(r, c+1)))
                LL = min(J(L(r,c)), J(L(r, c+1))) ;

                J1 = J(L(r,c));
                J2 = J(L(r, c+1)) ;
                for k=1:N
                    if (J(k) == J1 | J(k) == J2)
                        J(k) = LL ;
                    end ;
                end;

            end ;
            if (L(r+1, c) ~= 0 & J(L(r,c)) ~= J(L(r+1, c)))
                LL = min(J(L(r,c)), J(L(r+1, c))) ;
                J1 = J(L(r,c)) ;
                J2 = J(L(r+1, c)) ;
                for k=1:N
                    if (J(k) == J1 | J(k) == J2)
                        J(k) = LL ;
                    end ;
                end;
            end ;
        end ;
    end ;
end ;

for r = 1:S(1),
    for c = 1:S(2),
        if (L(r,c) ~= 0)
            L(r,c) = J(L(r,c)) ;
        end ;
    end ;
end ;
  
%      figure(2)
%      imagesc(L')
%      colormap('gray')
%      axis image
%      pause(0.1)

