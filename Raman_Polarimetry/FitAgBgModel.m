function result = FitAgBgModel(theta,parallel,cross,settings)
%% FitAgBgModel.m fits Raman data taking into account both parallel-polarized
%% cross-polarized channels.

    if nargin<4, settings = struct; end
    defaults = struct('ratioBounds',[-20,20],'ratioGridPoints',2001, ...
        'optimOptions',optimset(Display="off",TolX=1e-10,MaxIter=1000,MaxFunEvals=3000));
    for name = string(fieldnames(defaults)).'
        if ~isfield(settings,name), settings.(name) = defaults.(name); end
    end
    [nT,~,nModes] = size(parallel);
    parallel = double(parallel);
    cross = double(cross);
    theta = double(theta);
    data = {parallel,cross};
    traces = cell(nT,nModes,2);
    for g = 1:2
        for t = 1:nT
            for m = 1:nModes
                y = reshape(data{g}(t,:,m),[],1);
                valid = isfinite(y);
                y = y(valid);
                traces{t,m,g} = struct('y',y,'valid',valid,'sst',sum((y-mean(y)).^2));
            end
        end
    end
    cos2 = cos(theta).^2;
    sin2 = sin(theta).^2;
    ratio = zeros(1,nModes);
    for m = 1:nModes
        objective = @(r) shapeLoss(r,cos2,sin2,traces(:,m,1));
        ratio(m) = searchRatio(objective,settings);
    end
    coefAg = zeros(nT,nModes,2);
    coefBg = coefAg;
    predicted = {nan(size(parallel)),nan(size(cross))};
    for m = 1:nModes
        for g = 1:2
            for t = 1:nT
                trace = traces{t,m,g};
                design = tensorBasis(cos2(:,t),sin2(:,t),ratio(m),g);
                coef = nonnegativeFit(design(trace.valid,:),trace.y);
                predicted{g}(t,:,m) = (design*coef).';
                coefAg(t,m,g) = coef(1);
                coefBg(t,m,g) = coef(2);
            end
        end
    end
    fractionAg = coefAg./(coefAg+coefBg);
    fractionAg(coefAg+coefBg==0) = 0.5;
    result.ratio = ratio;
    suffix = ["Parallel","Cross"];
    for g = 1:2
        result.("coefAg"+suffix(g)) = coefAg(:,:,g);
        result.("coefBg"+suffix(g)) = coefBg(:,:,g);
        result.("weightAg"+suffix(g)) = fractionAg(:,:,g);
        result.("predicted"+suffix(g)) = predicted{g};
    end
end

function loss = shapeLoss(r,cos2,sin2,traces)
    loss = 0;
    for t = 1:numel(traces)
        trace = traces{t};
        design = tensorBasis(cos2(trace.valid,t),sin2(trace.valid,t),r,1);
        residual = trace.y-design*nonnegativeFit(design,trace.y);
        loss = loss+sum(residual.^2)/trace.sst;
    end
end

function ratio = searchRatio(objective,settings)
    bounds = settings.ratioBounds;
    grid = linspace(bounds(1),bounds(2),settings.ratioGridPoints);
    if bounds(1)<0 && bounds(2)>0, grid = unique([grid,0]); end
    gridLoss = arrayfun(objective,grid);
    localMinimum = find(gridLoss(2:end-1)<=gridLoss(1:end-2) & ...
        gridLoss(2:end-1)<=gridLoss(3:end))+1;
    count = 2+numel(localMinimum);
    candidates = zeros(1,count+2*numel(localMinimum));
    losses = inf(size(candidates));
    candidates(1:count) = [grid(1),grid(end),grid(localMinimum)];
    losses(1:count) = [gridLoss(1),gridLoss(end),gridLoss(localMinimum)];
    flat = max(gridLoss)-min(gridLoss)<=1e-12*max(1,max(abs(gridLoss)));
    if ~flat && ~isempty(localMinimum)
        for index = localMinimum([true,diff(localMinimum)>1])
            edges = unique([grid(index-1),max(grid(index-1),min(0,grid(index+1))),grid(index+1)]);
            for j = 1:numel(edges)-1
                count = count+1;
                [candidates(count),losses(count)] = fminbnd(objective,edges(j),edges(j+1),settings.optimOptions);
            end
        end
    end
    [~,best] = min(losses(1:count));
    ratio = candidates(best);
end

function design = tensorBasis(cos2,sin2,r,geometry)
    if geometry==1
        design = [(cos2+r*sin2).^2,4*cos2.*sin2];
    else
        design = [(r-1)^2*cos2.*sin2,(cos2-sin2).^2];
    end
end

function coefficients = nonnegativeFit(design,observed)
    scales = vecnorm(design);
    scales(scales==0) = 1;
    normalized = design./scales;
    candidates = [zeros(2,1),diag(max(normalized.'*observed,0))];
    unconstrained = pinv(normalized)*observed;
    if all(unconstrained>=0), candidates(:,1) = unconstrained; end
    [~,best] = min(sum((normalized*candidates-observed).^2,1));
    coefficients = candidates(:,best)./scales.';
end
