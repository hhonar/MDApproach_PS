function [rho] = circularslidingwindow(alpha, beta, w,method)
rho = NaN(size(alpha));

    switch method
        case {'vonmises'}
            kappa = pi/16;%
            Win = vonmissespdf(linspace(-pi,pi,w),0,kappa,1)';
            Win = Win/sum(Win);
        otherwise
            %disp('not tapered');
            Win = 1;
    end

    for i = 1:length(alpha)-w
        rho(i+ceil(w/2)) = circ_corrcc(alpha(i:i+w-1).*Win,beta(i:i+w-1).*Win);
    end

end
