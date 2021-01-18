function details = Cont_alpha_1(hFun,u0,lambda0,LimSteps,varargin)
% REWRITING CONTINUATIN
% 
% November 2018 -- implemented MoveMesh, which takes in u and lambda, 
%                   changes them and returns them
% NOTE: 
% To set any global vars inside hFun to those corresponding
% to accepted u one must call it either in AssingOutputs here or
% call hFun from inside the monitor function inside the function calling 
% contHDFbeta1. This is because of the computation of jacobian: the last call to 
% hfun is with lambda+eps, which changes any global vars depending on
% lambda to those corresponding to lambda+eps
% 
% TODO (* -- absolutely asap, becasue this slows down)
% *. Remove dynamically growing arrays, such as details. Return number of
%    steps instead
% 
% *. Pass step index into monitor function as a third argument. This can
%     be used with global arrays to record data.
% 
% 1. FIX RE-EVALUATING EXPRESSIONS USED IN hFun IN MonitorFun
%    E.g., allow varargout in hFun and varargin in MonitorFun to pass
%    these. Modify FunJac and FunJacAna accordingly. Modify also calls to
%    them. 
%    Now using global variables as a work-around. Consider persistent vars.
% 
% 2. A lot of solution measures are being computed inside MonitorFun. Find
%    a way to cumumlatively save this data into variables during run times.
%    For 2D computations, this should be dumped on hard drive after every
%    N0 calls. Use matfile object
% 
% 3. Remove .hdf file use and use .mat files instead
% 
% 4. When BC is contact with reservoir, we need to change amp to the new
%    value of rhoB. Doing it inside target fun is very costly, so now I am
%    doing it inside monitorfun, essentially using the bath from the
%    previous mu. Think how this can be automated. Maybe, implement PRE-PROCESSOR 
%    function to set vals of some vars, on which target fun depends, and 
%    then call target fun from inside the pre-processor fun.
% 
% 5. No need to return anything from monitor fun. Turns out the best practice
%    is to dump into file from monitor fun. Maybe a better practice is to 
%    return stuff from monitor fun via varargout, which goes inside POST-PROCESSOR,
%    where stuff gets dumped into file. This will de-mess monitor fun.
% 
% 6. JUMPED TOO FAR OPTION DOES NOT WORK, PROBABLY, 
%    FOR THE SAME REASON AS STOP BUTTON
%    STOPS ON NEXT ITERATION, NOT ON CURRENT.
% 
% 7. Figure out, how I can change continuation parameters and use the same
%    MonitorFcn. -- See rain_model_D
% 
% ULTIMATE VERSION OF contHDF, previous is contHDFbeta
% NOTE: TEST ADAPTIVE WITH BOTH, KELLER AND MOORE-PENROSE
% MODES: 
%        start        - get first (N+1)-point by contHDF engine, treat u0,
%                       lambda0 as IG
% 
%        proceed      - first (N+1)-point supplied, get tangent at it and
%                       proceed, treat u0, lamda0 as 1st (N+1)-point
% 
%        appendhdf    - continue with previous continuation, disregard
%                       input u0, lambda0, read those and tangent from last 
%                       dataset of file
% to implement:
%        revappendhdf - continue with previous continuation, but inverting
%                       the direction of change of lambda0, starting from
%                       the first dataset
    % ARGUMENT CHECKS
    if mod(nargin,2)>0, error('Mismatched number of input arguments'); end

    lastwarn('');
    epsilon = sqrt(eps);
    
    % DEFAULT OPTIONS
    opts = struct('movemeshfcn',[],...
                  'monitorfcn',[],...
                  'step',0.1,...
                  'adaptive',[],...
                  'maxstep',inf,...
                  'minstep',eps,...
                  'jacobian',[],...
                  'algorithm','moore-penrose',...
                  'maxiter',4*length(u0),...
                  'direction',1,...
                  'tolfun',1e-8,...
                  'filename','',...
                  'mode','start',...
                  'stability','off',...
                  'names',{{'U','lambda'}},...
                  'display','on',...
                  'stopbtn',[]);
    
    details = []; MonVal = [];
              
    % READ OPTIONS
    for i=1:2:nargin-4,
        opts.(lower(varargin{i}))=varargin{i+1};
    end
        
    
    if strcmpi(opts.display,'on'), isDisp = 1;
    else isDisp = 0;
    end
              
    if ~isempty(opts.stopbtn), hStopBtn = []; end 
    
    ds = opts.step;   
    hMonitor = opts.monitorfcn;
hMoveMesh = opts.movemeshfcn;    
    TolFun = opts.tolfun;
    maxiter = opts.maxiter;
    N = length(u0);
    
    if ~isempty(hMonitor), isMon = 1; else isMon = 0; end
if ~isempty(hMoveMesh), isMoveMesh = 1; else isMoveMesh = 0; end
        
    if strcmpi(opts.stability,'on'), isStab = 1; else isStab = 0; end;
    
    % CHECK ADAPTIVITY
    if isempty(opts.adaptive)
        isAdaptive = 0;
    else
        isAdaptive = 1;
        kap0    = opts.adaptive(1);
        al0     = opts.adaptive(2);
        delt0    = opts.adaptive(3);
    end
      
    
    % CHECK JACOBIAN INPUT
    if strcmpi(opts.jacobian,'on'), hFunJac = @FunJacAna;
    else                            hFunJac = @FunJac;
    end
    
    % u0, lambda0 treated as IG    
    if strcmpi(opts.mode,'start')
        % SOLVE FOR THE FIRST TIME & GET JACOBIAN
        [f,Jac] = hFunJac(u0,lambda0);
        isOK = 0;
        
        for i = 1:maxiter,
            % CHECK CONVERGENCE
            if(norm(f,2) <= TolFun), isOK = 1; break; end;
        
            % NEWTON STEP
            u0 = u0 - (Jac(:,1:N)\f);
                
            % NEW FUNCTION & JACOBIAN VALUES
            [f,Jac] = hFunJac(u0,lambda0);
        end
    
        if isOK, reject = 0; iStep = 1;
        else error('MATLAB:CONT','Bad initial guess'); 
        end
        
        u0 = [u0;lambda0]; u = u0; 
        if isAdaptive, kap = kap0; delt = delt0; al = al0; end
        
        
    
    % CONTINUE A ONCE STARTED CONTINUATION, REQUIRES FILE opts.filename.hdf
    elseif  strcmpi(opts.mode,'appendhdf') 
        iStep = 0; reject = 0;
        hinfo = hdf5info(opts.filename);
        DataSets = size(hinfo.GroupHierarchy.Groups,2)-1;
        u0 = hdf5read(opts.filename,[num2str(DataSets),'/',opts.names{1}]);
        lambda0 = hdf5read(opts.filename,[num2str(DataSets),'/',opts.names{2}]);
        t0 = hdf5read(opts.filename,[num2str(DataSets),'/tau']);
        if opts.direction == -1,  t0 = -t0;  end
        u0 = [u0;lambda0]; u = u0; 
        if isAdaptive, kap = kap0; delt = delt0; al = al0; end
    
    % CONTINUE A KNOWN SOLUTION, u0, lambda0 treated as solution, not as IG    
    elseif strcmpi(opts.mode,'proceed')
        iStep = 0; reject = 0;
        [f,Jac] = hFunJac(u0,lambda0); %#ok
        u0 = [u0;lambda0]; u = u0;
    else error('ContHdf: unknown mode')
    end
   
    if ~strcmpi(opts.mode,'appendhdf')
        % GET TANGENT AS ZERO VECTOR OF Jac
        % add "some" row at the end of Jac, to make it square and solve by "\"
        t0 = [Jac;ones(1,N+1)]\[zeros(N,1); 1;];
        t0 = t0/norm(t0,2);
        % if "some" row produced a zero det, do it proper
        if sum(isnan(t0))>0 || sum(isinf(t0))>0,
            fprintf('contHDF warning: computing initial dz by calling null()\n')
            keyboard
            t0 = null(Jac); % this should not happen. Maybe, the continuatin parameter chosen is part of the solution, and not a free parameter
            
        end
        % chose tangent projection on parameter axis with sign given by 
        % opts.direction 
        t0 = (opts.direction*sign(t0(end)+eps))*t0;

        t = t0; 
        AssignOutputs; 
   end
    
    % CONTINUE VIA Penrose-Moore or Keller
    switch lower(opts.algorithm),
        case 'moore-penrose'
            while iStep<LimSteps,
                
                % restart timer
                if isDisp, tic; end
                
                % count converged steps
                if ~reject, iStep = iStep + 1; end;
                
                % INITIALIZE FLAGS
                isTerminal = 0; % 1 means terminate continuation
                reject = 1;     % 1 means re-do Newton from u0
                isOK = 0;       % 1 means Newton converged
                lastwarn('');   % clear last warning
                
                % PREDICTOR STEP
                t = t0;
                u = u0 + t0*ds;
       
                % NEWTON ITERATION
                for NewtonStep = 1:maxiter,
                    
                    % EVALUATE F & JACOBIAN
                    [f,Jac] = hFunJac(u(1:end-1),u(end));

                    % EXIT CURRENT NEWTON BY CONVERGENCE
%                     if(norm(f,2) <= TolFun), isOK = 1; break; end
                    if(norm(f,2) <= TolFun)&&NewtonStep>1, isOK = 1; break; end % when converged on 1st iteration, al,delt,kap were undefined
                    
                    % EXIT CURRENT NEWTON BY IDENTIFYING IT AS 
                    % NON-CONVERGENT:
                    % 1. something wrong: we entered a complex plane
                    % 2. MATLAB has issued a warning
                    if norm(imag(f),2)>TolFun || ...
                       ~isempty(lastwarn), break; end
                    
                    % CORRECTOR
                    delta = -[Jac ; t']\[f Jac*t; 0 0];
                    u = u + delta(:,1);
                    t = t + delta(:,2);
                    t = t/norm(t,2);  
                    
if NewtonStep == 1 && isAdaptive,
    Jac1 = Jac;
    al = acos(dot(t,t0));
    delt = norm(delta(:,1),2); 
    kap = kap/2; %increase kap in case Newton converges in 1 iteration
end
if NewtonStep ==2 && isAdaptive,
    kap = norm(Jac1\f,2)/delt;
    clear Jac1
%     fprintf(1,'kap = %6.4f, al = %6.4f, delt = %6.4f \n',kap,al,delt);
end
                    clear('Jac');
                 end
   
                % POST PROCESSING: 1 -- finish cont
                if PostProcessing, break; end;
            end
        case 'keller'
             error('Not emplemeted in new version yet')
        otherwise
            error('MATLAB:CONT',['Invalid algorithm.',...
                ' Only ''Moore-Penrose'' and ''Keller'' are acceptable']);
    end
 


          
    % ASSIGNS OUTPUTS
    function AssignOutputs
        
    % If using global vars inside hFun, their vals depending on lambda, [equivalently, u(end)] 
    % would be correspodning to u(end)+eps due to computation of the jacobian in FunJacAna. 
    % To set global vars correct, either call hFun inside the monitor function or uncomment the line below:
    % hFun(u(1:end-1),u(end));  
    

% here change u between successful iterations    
if isMoveMesh,
    [u1,lam1] = hMoveMesh(u(1:end-1),u(end));
    u = [u1;lam1];
end
    
        if isMon,
  
            [MonVal,status] = hMonitor(u(1:end-1),u(end));
            
            % CONTROL CONTINUATION FROM INSIDE MONITOR FCN
            % terminate (e.g. enough solutions obtained)
            if isnan(status),isTerminal = 1; end
             % reject (e.g., accidentally jumped too far in measure)
            if isinf(status), reject = 1; end
            
            % ---------------------------------
            % STOP BUTTON 
            % ---------------------------------
            if ~isempty(opts.stopbtn)

                % identify figure to put button on:
                % 0     --  assume there is only one figure and use it, 
                % else  --  assume that figure handle is passed
                if opts.stopbtn==0, hFig = findobj('type','figure');
                else hFig = opts.stopbtn;
                end

                if isempty(hStopBtn) % must assign only once, or handle changes
                   % delete previous button ?
                   delete(findobj('tag','ContStop'));
                   % position button on figure
                   pos = get(0,'DefaultUicontrolPosition');
                   pos(1) = pos(1) - 15;            
                   pos(2) = pos(2) - 15;    
                   % form button
                   hStopBtn = uicontrol('Parent',hFig,...
                                        'Style','pushbutton', ...         
                                        'Callback',@StopBtnCallback,...
                                        'String','Stop', ...
                                        'Tag','ContStop',...
                                        'Position',pos);
                   set(hFig,'UserData',0);
                 end

                 % if pressed, terminate
                 if get(hStopBtn,'UserData'), isTerminal = 1; end
            end
        end
        
        % ----------------
        % WRITE TO FILE        
        % ----------------
        if ~reject &&...   % ensure file entries are "in order"
           ~isempty(opts.filename), % file output requested
    
            try
                hinfo = hdf5info(opts.filename);
                DataSets = ['/',num2str(size(hinfo.GroupHierarchy.Groups,2)),'/'];
                Md = 'append';
            catch %#ok
                DataSets = '/0/';
                Md = 'overwrite';
            end
           
            hdf5write(opts.filename,...
                      [DataSets,opts.names{1}],u(1:end-1),...
                      [DataSets,opts.names{2}],u(end),...
                      [DataSets,'tau'],t,...
                      'WriteMode',Md);
                  
            % LINEAR STABILITY               
            if isStab, hdf5write(opts.filename,...
                                 [DataSets,'stability'],...
                                 isStable(Jac),'WriteMode','append');
            end
            
            % MONITOR FUNCTION OUTPUT
            if isMon, hdf5write(opts.filename,...
                                [DataSets,'monitor'],MonVal,...
                                'WriteMode','append');
            end
        end        
    end

    % RETURNS TARGET FUNCTION AND NUMERICAL JACOBIAN
    function [FF,JJ] = FunJac(x,lambda)
       FF = hFun(x,lambda);
       
       % NUMERICAL JACOBIAN
       x = [x;lambda];
       JJ = zeros(N,N+1);
       for n=1:N+1,
           x_new = x;  x_new(n) = x(n) + epsilon;
           f_new = hFun(x_new(1:N),x_new(N+1));
           JJ(:,n) = (f_new-FF)/epsilon;
       end
    end

    % RETURNS TARGET FUNCTION AND USER-SUPPLIED JACOBIAN
    function [FF,JJ] = FunJacAna(x,lambda)
        [FF,JJ] = hFun(x,lambda);
        
        
        % this step changes values of global vars inside hfun depending on
        % lambda, so there value does not correspond to accepted u, unless
        % there is an additional call to hfun in AssignOutputs
        if size(JJ,2)==N,
           JJ = [JJ, (hFun(x,lambda + epsilon) - FF)/epsilon]; 
        end
    end

    % POST PROCESSING (now value of reject is 1)
    function out = PostProcessing
        out = 0;  % do we finish with continuation ? 
        
        if isAdaptive || isDisp
            alpha = acos(dot(t,t0)); % cos of angle between tangents at u, u0
            dX = norm(u-u0,2); % distance between solutions,? /norm(u0) ?
        end
        
        % -----------------------
        %   ADAPTIVE STEP
        % -----------------------
        if isAdaptive,
            
            % can not further reduce step
            if ds < opts.minstep, 
                isTerminal = 1;
                warning('MATLAB:CONT','Smallest step size reached');
            end
 
%             % CHANGE STEP
%             
%             % exceeded tolAlph
%             % reject and reduce step based on angle 
%             % (turn of curve ?, special for alpha <0?)
%             if alpha<tolAlph
%                     ds = ds*min((1-abs(alpha))^2,...
%                                abs(alpha)^2); 
%                     if isDisp, fprintf(1,'tolAlph\n'); end
%             % exceeded tolX
%             % reject and half step
%             elseif dX>tolX
%                     ds = ds/2; 
%                     if isDisp, fprintf(1,'tolX\n'); end
%             % converged, speed-up               
%             elseif isOK && NewtonStep<Nthresh,
%                    ds = min(ds*(1+alpha),opts.maxstep); 
%                    reject = 0;
%             
%             % non-converged
%             % reject and reduce step by at least half
%             elseif ~isOK,
%                    ds = ds*min([(1-abs(alpha))^2,...
%                                abs(alpha)^2,...
%                                .5]); 
%                    if isDisp, fprintf(1,'~isOK\n'); end
%             % Newton converged, but at more than Nthresh steps
%             else
%                 reject = 0;
%             end
 
reject = 0;
% kap0  = .8;
% al0   = .1;
% delt0 = .1;
[sqrt(kap/kap0),sqrt(delt/delt0),al/al0]
cc = max([sqrt(kap/kap0),sqrt(delt/delt0),al/al0]);
cc = max([min([cc,2]),.5]); % .5 <cc< 2 
ds = min([opts.maxstep,ds/cc]);

if cc==2 , reject = 1; end



        % ----------------------
        %   NON-ADAPTIVE
        % ----------------------
        % converged
        elseif isOK, reject = 0;
        
        % non-converged
        else isTerminal = 1;
            
        end
        
        % -------------------------------
        %   PROCESS reject & isTerminal
        % -------------------------------
        % store if not rejected, 
        % exit by output from monitorfcn -- changes isTerminal
        % possibility of rejection by jumping too far -- changes reject
        if ~reject, AssignOutputs; end
        
        if ~reject, u0 = u; t0 = t; 
%         else ds = ds/2; % jumped too far in measure
        end
      
        % ---------------------------------------
        %   CONSOLE OUTPUT
        % ---------------------------------------
        if isDisp && isAdaptive && ~isempty(opts.filename),
            details = [details;...
                   [iStep,...
                   ~reject,...
                   NewtonStep,...
                   dX,...
                   alpha,...
                   ds,...
                   kap,...
                   delt,...
                   al,...
                   1/cc,...
                   toc]]; 
            if mod(iStep,10)==0 || iStep <= 2,
                fprintf(1,'|---------------------------------------------------------------|\n');
                fprintf(1,'| # |Stps|  dX  |alpha |  ds  |kappa |delta |alpha0|koeff |time |\n');
                fprintf(1,'|---------------------------------------------------------------|\n');
            end
            fprintf(1,'|%3i|%1i/%2i|%6.4f|%6.4f|%6.4f|%6.4f|%6.4f|%6.4f|%6.4f|%5.3f|\n',details(end,:))
        
        elseif isDisp
            
                if mod(iStep,10)==0 || iStep  <= 2,
                fprintf(1,'|---------------------------------------------------------------|\n');
                fprintf(1,'| # |Stps|  dX  |alpha |  ds  |time |\n');
                fprintf(1,'|---------------------------------------------------------------|\n');
                end
                fprintf(1,'|%3i|%1i/%2i|%6.4f|%6.4f|%6.4f|%5.3f|\n',...
                        [iStep,~reject,NewtonStep,dX,alpha,ds,toc]); 
        end
        
        if isempty(opts.filename)
            details = MonVal;
        end
       
        
          
        if isTerminal,
            iStep = iStep-(~isOK); % discard index of non-converged step, 
                                   % which caused temination
            warning('MATLAB:CONT','Continuation terminated.');
            out = 1;
        end
    end
end

% IF STOP BUTTON PRESSED
function StopBtnCallback(gco,varargin)
    % modify property of button object, which is later interpreted
    set(gco,'UserData',1);
end

% DETERMINES STABILITY
function out = isStable(Jac)
    eigvals = real(eig(Jac));
    if all(eigvals>0), out = 1; 
    elseif all(eigvals<0), out = -1;
    else out = 0; 
    end
end