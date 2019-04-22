function [z,m,s] = nanzscore(x,opt,varargin)
    if (nargin < 3) % check input
        dim = find(size(x)>1,1);
        if isempty(dim)
            dim = 1;
        end
    else
        dim = varargin{1};
    end

    if (nargin < 2) || isempty(opt)
        opt = 0;
    end

    % compute mean value(s) and standard deviation(s)
    m = nanmean(x,dim);
    s = nanstd(x,opt,dim);    
    % computer z scores
    z = (x-repmat(m,size(x)./size(m)))./repmat(s,size(x)./size(s));
end


function m = nanmean(x, varargin)
    if (nargin < 2) % check input
        dim = find(size(x)>1,1);
        if isempty(dim)
            dim = 1;
        end
    else
        dim = varargin{1};
    end

    % determine number of regular (not nan) data points
    n = sum(~isnan(x),dim);

    % replace nans with zeros and compute mean value(s)
    x(isnan(x)) = 0;
    n(n==0) = nan;
    m = sum(x,dim)./n;
end


function s = nanstd(x, opt, varargin)
    if (nargin < 3) % check input
        dim = find(size(x)>1,1);
        if isempty(dim)
            dim = 1;
        end
    else
        dim = varargin{1};
    end
     
    if (nargin < 2) || isempty(opt)
        opt = 0;
    end
 
    % determine number of regular (not nan) data points and nans
    n = sum(~isnan(x),dim);
	nnan = sum(isnan(x),dim);
     
    % replace nans with zeros, remove mean value(s) and compute squared sums
    x(isnan(x)) = 0;
    m = sum(x, dim)./n;
    x = x-repmat(m, size(x)./size(m));
    s = sum(x.^2, dim);

    % remove contributions of added zeros
    s = s-(m.^2).*nnan;

    % normalization
    if (opt == 0)
        s = sqrt(s./max(n-1,1));
    elseif (opt == 1)
        s = sqrt(s./n);
    else
        error('unkown normalization type');
    end
end