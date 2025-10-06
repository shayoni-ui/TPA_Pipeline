classdef NCA
   %NCA a controlled object for non-compartmental analysis
   
   properties
       Regimen      = '';
       RegimenName  = '';
       Analyte      = '';
       Species      = '';
       Group        = '';
       Sex          = '';
       Subject      = '';
       Period       = '';
       Tissue       = '';
       Tmax         = NaN;
       TrueTmax     = false;
       Cmax         = NaN;
       Tlast        = NaN;
       Clast        = NaN;
       Clast_pred   = NaN;
       AUC          = NaN;
       AUCinf       = NaN;
       AUCinf_pred  = NaN;
       AUMC         = NaN;
       AUMCinf      = NaN;
       AUMCinf_pred = NaN;
       MRT          = NaN;
       MRTinf       = NaN;
       MRTinf_pred  = NaN;
       lambdaZ      = NaN;
       lambdaZ_pts  = NaN;
       lambdaZ_rng  = [NaN NaN];
       lambdaZ_nHL  = NaN;  % Number of half-lives used
       % to calculate lambdaZ
       Rsq          = NaN;
       Rsq_adjusted = NaN;
       Vss          = NaN;
       Vss_pred     = NaN;
       Vz           = NaN;
       Vz_pred      = NaN;
       CL           = NaN;
       CL_pred      = NaN;
       t12          = NaN;
       t12_CI       = [NaN NaN];
       Kp           = NaN;
       Kp_method    = '';
       Kb           = NaN;
       Kb_method    = '';
       LLQ          = 0;
       LLQmethod    = 'ignore';
       Dose         = NaN;
       DoseUnit     = 'mg';
       DoseLength   = 0;
       Data         = struct('time',[],'conc',[],'mdl',[]);
       AmountUnit   = 'ng';
       VolumeUnit   = 'mL';
       TimeUnit     = 'hr';
       Comments     = {};
       
   end
   
   % Class methods
   methods
       function obj = NCA(nva)
          arguments
               nva.Regimen      = '';
               nva.RegimenName  = '';
               nva.Analyte      = '';
               nva.Species      = '';
               nva.Group        = '';
               nva.Sex          = '';
               nva.Subject      = '';
               nva.Period       = '';
               nva.Tissue       = '';
               nva.Tmax         = NaN;
               nva.TrueTmax     = false;
               nva.Cmax         = NaN;
               nva.Tlast        = NaN;
               nva.Clast        = NaN;
               nva.Clast_pred   = NaN;
               nva.AUC          = NaN;
               nva.AUCinf       = NaN;
               nva.AUCinf_pred  = NaN;
               nva.AUMC         = NaN;
               nva.AUMCinf      = NaN;
               nva.AUMCinf_pred = NaN;
               nva.MRT          = NaN;
               nva.MRTinf       = NaN;
               nva.MRTinf_pred  = NaN;
               nva.lambdaZ      = NaN;
               nva.lambdaZ_pts  = NaN;
               nva.lambdaZ_rng  = [NaN NaN];
               nva.lambdaZ_nHL  = NaN;  % Number of half-lives used
                   % to calculate lambdaZ
               nva.Rsq          = NaN;
               nva.Rsq_adjusted = NaN;
               nva.Vss          = NaN;
               nva.Vss_pred     = NaN;
               nva.Vz           = NaN;
               nva.Vz_pred      = NaN;
               nva.CL           = NaN;
               nva.CL_pred      = NaN;
               nva.t12          = NaN;
               nva.t12_CI       = [NaN NaN];
               nva.Kp           = NaN;
               nva.Kp_method    = '';
               nva.Kb           = NaN;
               nva.Kb_method    = '';
               nva.LLQ          = 0;
               nva.LLQmethod    = 'ignore';
               nva.Dose         = NaN;
               nva.DoseUnit     = 'mg';
               nva.DoseLength   = 0;
               nva.Data         = struct('time',[],'conc',[],'mdl',[]);
               nva.AmountUnit   = 'ng';
               nva.VolumeUnit   = 'mL';
               nva.TimeUnit     = 'hr';
               nva.Comments     = {};
          end
          flds = fields(obj);
          for ii = 1:length(fields(nva))
            fld = fields(nva);
            if ~contains(flds, fld(ii))
                error('Unrecognized input "%s"', fld);
            else
                obj.(fld{ii}) = nva.(fld{ii});
            end
          end
%           for i = 1:2:length(varargin)
%               if sum(strcmpi(flds,varargin{i}))
%                   fld       = flds{strcmpi(flds,varargin{i})};
%                   obj.(fld) = typecheck(varargin{i+1}, ...
%                       class(obj.(fld)),varargin{i});
%               elseif endsWith(varargin{i},'s') &&...
%                       sum(strcmpi(flds,varargin{i}(1:length(varargin{i}) ...
%                       -1))) % Replace trailing "s" and try again
%                   fld       = flds{strcmpi(flds, ...
%                       varargin{i}(1:length(varargin{i})-1))};
%                   obj.(fld) = typecheck(varargin{i+1}, ...
%                       class(obj.(fld)),varargin{i});
%               else, error('Unrecognized input "%s"',varargin{i})
%               end
%           end
      end
      function disp(obj)
            s=jsonencode(obj,'PrettyPrint',true);
            fprintf('%s\n',s);
%           regimens = unique({obj.Regimen});
% 
%           for i = 1:length(regimens)
%               dispobj = obj(strcmpi({obj.Regimen},regimens{i}));
%               amt_unit = unique({dispobj.AmountUnit});
%               vol_unit = unique({dispobj.VolumeUnit});
%               tim_unit = unique({obj.TimeUnit});
%               concunit = sprintf('%s/%s',amt_unit{:},vol_unit{:});
%               aucunit  = sprintf('%s*%s/%s',amt_unit{:}, ...
%                   tim_unit{:},vol_unit{:});
%               
%               fprintf('NCA for %s\n',regimens{i});
%               fprintf(['\t%-10.10s\n\t%-5.5s\n\t%-10.10s\n\t%-12.12s\n' ...
%                   '\t%-12.12s\n\t %-7.7s\n\t %-10.10s\n\t%-10.10s\n' ...
%                   '\t %-11.11s\n'],'', ...
%                   '','','','','Cmax','AUC_last','AUC_inf','Half-Life');
%               fprintf(['\t%-10.10s\n\t%-5.5s\n\t%-10.10s\n\t%-12.12s\n' ...
%                   '\t%-12.12s\n\t %-7.7s\n\t %-10.10s\n\t%-10.10s\n' ...
%                   '\t%-8.8s\n'], ...
%                   'Subject','Sex','Analyte','Period','Tissue', ...
%                   concunit,aucunit,aucunit,tim_unit{:});
%               for j = 1:length(dispobj)
%                   fprintf('\t%-10.10s\n',   dispobj(j).Subject);
%                   fprintf('\t%-5.5s\n',     dispobj(j).Sex);
%                   fprintf('\t%-10.10s\n',   dispobj(j).Analyte);
%                   fprintf('\t%-12.12s\n',   dispobj(j).Period);
%                   fprintf('\t%-12.12s\n',   dispobj(j).Tissue);
%                   fprintf('\t%8.3g\n',      dispobj(j).Cmax);
%                   fprintf('\t %8.3g\n',     dispobj(j).AUC);
%                   fprintf('\t%8.3g\n',      dispobj(j).AUCinf);
%                   fprintf('\t   %3.1f\n',   dispobj(j).t12);
%                   fprintf('\n');
%               end
%           end
      end
      function strc = nca2struct(obj)
          flds = fields(obj);
          for i = 1:numel(obj)
            for j = 1:length(flds)
                strc(i).(flds{j}) = obj(i).(flds{j}); %#ok<*AGROW>
            end
          end
      end
      
      function plot(obj,ax)
          if length(obj) > 1
              if nargin ==1, figure; end
              for i = 1:length(obj)
                  if nargin ==1
                      ax = subplot(round(sqrt(length(obj))), ...
                          ceil(sqrt(length(obj))),i); 
                  else
                      set(ax,'yscale','log')
                  end
                  plot(obj(i),ax);
              end
              return
          end
          
          if nargin == 1, figure; ax = subplot(1,1,1); end
          
          if ~isempty(obj.Data.mdl)
            t_mdl = linspace(0,max(obj.Data.time),500);
            c_mdl = exp(obj.Data.mdl.Coefficients.Estimate(1) + ...
                obj.Data.mdl.Coefficients.Estimate(2) .* t_mdl);
          else
            t_mdl = [];
            c_mdl = [];
          end
          
          if ~isempty(obj.Data.time) && ~isnan(obj.lambdaZ_pts)
              t_fit = obj.Data.time(end-obj.lambdaZ_pts+1:end);
              c_fit = obj.Data.conc(end-obj.lambdaZ_pts+1:end);

              t_oth = obj.Data.time(1:end-obj.lambdaZ_pts);
              c_oth = obj.Data.conc(1:end-obj.lambdaZ_pts);
          else
              t_fit = [];
              c_fit = [];
              t_oth = obj.Data.time;
              c_oth = obj.Data.conc;
          end
          
          plot(ax,t_oth,c_oth,'bx',t_fit,c_fit,'rx',t_mdl,c_mdl,'r-')
          set(ax,'yscale','log')
          xlabel(sprintf('Time, %s',obj.TimeUnit))
          ylabel(sprintf('Conc, %s/%s',obj.AmountUnit,obj.VolumeUnit))
          switch lower(obj.TimeUnit)
              case {'hr','hour','hours'}
                  if max(obj.Data.time) <= 6            
                      ticks = 0:max(obj.Data.time);
                  elseif max(obj.Data.time) <= 36       
                      ticks = 0:4:max(obj.Data.time);
                  elseif max(obj.Data.time) <= 96       
                      ticks = 0:12:max(obj.Data.time);
                  elseif max(obj.Data.time) <= 168      
                      ticks = 0:24:max(obj.Data.time);
                  elseif max(obj.Data.time) <= 336      
                      ticks = 0:48:max(obj.Data.time);
                  elseif max(obj.Data.time) <= 24*28*2 
                      ticks = 0:168:max(obj.Data.time);
                  elseif max(obj.Data.time) <= 24*28*4 
                      ticks = 0:336:max(obj.Data.time);
                  elseif max(obj.Data.time) <= 24*28*12
                      ticks = 0:672:max(obj.Data.time);
                  else, ticks = ax.XTick;
                  end
              case {'day','days'}
                  if max(obj.Data.time) < 7      
                      ticks = 0:max(obj.Data.time);
                  elseif max(obj.Data.time) < 60 
                      ticks = 0:7:max(obj.Data.time);
                  elseif max(obj.Data.time) < 365 
                      ticks = 0:28:max(obj.Data.time);
                  else, ticks = ax.XTick;
                  end
          end
          if isempty(ticks); ticks = ax.XTick; end
          ymin = floor(log10(ax.YLim(1)));
          ymax = ceil(log10(ax.YLim(2)));
          set(gca,'yscale','log','XTick',ticks, ...
              'YTick',10.^(ymin:ymax),'ylim',10.^[ymin ymax]);
          title({obj.Regimen; obj.Subject})
      end
      
    end 
end
