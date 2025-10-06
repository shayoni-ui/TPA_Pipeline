function [vnew, sp_input] = parameterizeDistribution(species, ...
    inputs, acatFlag)
% tag = ''; %For future parameterization of multi-drug models
if ~isempty(inputs.Tags)
    tag = ['_' inputs.Tags];
else
    tag = '';
end
%v = addModelVariant(species);
v = {};
%prmlimited = {};
phys       = inputs.Physiologies(strcmpi({inputs.Physiologies.Species}, ...
    species));
sp_input   = inputs.SpeciesSpecific(strcmpi({ ...
    inputs.SpeciesSpecific.Species},species));

vnew = sbiovariant(species);
vnew.addcontent({'parameter','wgt','value', ...
    phys.Properties.Weight});
vnew.addcontent({'parameter','mw','value',inputs.mw});
vnew.addcontent({'parameter',  ['B2P' tag], ...
        'Value',    sp_input.B2P});
if acatFlag
    vnew.addcontent({'parameter',  ...
            ['gut_fpe' tag],  'Value',    sp_input.gut_fpe});
    vnew.addcontent({'parameter', ['fu_ent' tag],'Value', ...
        sp_input.fu_ent});
end


%v_info = getModelValues(m);
%v = addVariantContent(v,{'parameter','wgt','value', ...
%    phys.Properties.Weight});
% v = addVariantContent(v,{'parameter','mw','value',inputs.mw});
% v = addVariantContent(v,{'Compartment','Depot',    'Capacity',  1});
% v = addVariantContent(v,{'parameter',  ['F' tag],        'Value',    1});
% v = addVariantContent(v,{'parameter',  ['ka' tag],       'Value',    1});

% prm = v_info(strcmpi(v_info(:,1),['t_diss' tag]),:);           
% if ~isempty(prm) 
%     v = addVariantContent(v,{'parameter', ['t_diss' tag],  ...
%         'Value', prm{2}}); 
% end
% prm = v_info(strcmpi(v_info(:,1),['Depot_Solubility' tag]),:); 
% if ~isempty(prm), v = addVariantContent(v,{'parameter', ...
%         ['Depot_Solubility' tag],  'Value', prm{2}}); end
% 
% prm = v_info(strcmpi(v_info(:,1),['B2P' tag]),:);     
% if ~isempty(prm), v = addVariantContent(v,{'parameter',  ['B2P' tag], ...
%         'Value',    sp_input.B2P});
% end
% prm = v_info(strcmpi(v_info(:,1),['gut_fpe' tag]),:);    
% if ~isempty(prm), v = addVariantContent(v,{'parameter',  ...
%         ['gut_fpe' tag],  'Value',    sp_input.gut_fpe}); 
% end

% prm = v_info(strcmpi(v_info(:,1),['fu_ent' tag]),:);            
% if ~isempty(prm), v = addVariantContent(v,{'parameter',  ['fu_ent' tag],  ...
%         'Value',    sp_input.fu_ent}); end


switch lower(inputs.fup_type)
    case {'adjusted','default'}   
%         v = addVariantContent(v,{'parameter',['fup' tag], ...
%             'Value',sp_input.adjusted_fup});
        vnew.addcontent({'parameter',['fup' tag], ...
            'Value',sp_input.adjusted_fup});
    case 'experimental'
%         v = addVariantContent(v,{'parameter', ...
%             ['fup' tag], 'Value',sp_input.fup});
        vnew.addcontent({'parameter',['fup' tag], 'Value',sp_input.fup});
    otherwise, error(['Invalid fup type, please use "Adjusted" ' ...
            'or "Experimental"'])
end

tis   = phys.Tissues;
switch lower(inputs.Distribution)
    case 'pbpk' % Add PBPK Distribiution and Clearance paramters
        % Assign tissue volumes
        for i = 1:length(tis)
            switch tis(i).ModelType
                case {'PerfusionLimited','PermeabilityLimited'}
                    vnew.addcontent({'compartment', ...
                            tis(i).Name,'Capacity',tis(i).Volume});

%                     cm = v_info(strcmpi(v_info(:,1),tis(i).Name),:);
%                     if ~isempty(cm)
%                         v = addVariantContent(v,{'compartment', ...
%                             tis(i).Name,'Capacity',tis(i).Volume});
%                         
%                     else, warning(['Unable to find "%s" in model,' ...
%                             ' compartmental and permeability limited PBPK ' ...
%                             'currently not supported'])
%                     end
                case 'BloodCompartment'
                    nm = strrep(strrep(tis(i).Name,' Supply',''), ...
                        ' Return','');
%                   v = addVariantContent(v,{'compartment',nm,'Capacity', ...
%                         tis(i).Volume});
                    vnew.addcontent({'compartment',nm,'Capacity', ...
                        tis(i).Volume});
                case 'GastroTissue' % No direct volume assigned to ACAT
                case 'PassThru'     % No direct volume assigned to 
                    % Hepatic Artery
                otherwise, error('Unrecongized tissue type "%s"', ...
                        tis(i).ModelType)
            end
        end
        
        % Add blood flow
%         v = addVariantContent(v,{'parameter','Q_Gut',          ...
%             'Value',tis(strcmp({tis.Name},'ACAT Gut')).TissuePerfusion});
        vnew.addcontent({'parameter','Q_Gut',          ...
            'Value',tis(strcmp({tis.Name},'ACAT Gut')).TissuePerfusion});
%         v = addVariantContent(v,{'parameter','Q_HepaticArtery', ...
%             'Value',tis(strcmp({tis.Name}, ...
%             'Hepatic Artery')).TissuePerfusion});
        vnew.addcontent({'parameter','Q_HepaticArtery', ...
            'Value',tis(strcmp({tis.Name}, ...
            'Hepatic Artery')).TissuePerfusion});
        tis = tis(strcmp({tis.ModelType}, ...
            'PerfusionLimited')|strcmp({tis.ModelType}, ...
            'PermeabilityLimited'));
        for i = 1:length(tis)
%             v = addVariantContent(v,{'parameter',['Q_' tis(i).Name], ...
%                 'Value',tis(i).TissuePerfusion});
            vnew.addcontent({'parameter',['Q_' tis(i).Name], ...
                'Value',tis(i).TissuePerfusion});
        end
        
        % Add Kp's and fut's
        kp_per = sp_input.Kp_perf;
        %kp_prm = sp_input.Kp_perm;
        for i = 1:length(tis)
            %cm = v_info(strcmpi(v_info(:,1),[tis(i).Name '_ec']),:);
            vnew.addcontent({'parameter',['Kp_' tis(i).Name tag],'Value', ...
                    kp_per.Kp(strcmpi({kp_per.Kp.Tissue},tis(i).Name)).Kp});
%             if isempty(cm)
%                 v = addVariantContent(v,{'parameter',[ ...
%                     'Kp_' tis(i).Name tag],'Value', ...
%                     kp_per.Kp(strcmpi({kp_per.Kp.Tissue},tis(i).Name)).Kp});
%             else
% 
%                 v = addVariantContent(v,{'parameter', ...
%                     ['Kp_' tis(i).Name tag],'' ...
%                     'Value',kp_prm.Kp(strcmpi({kp_prm.Kp.Tissue}, ...
%                     tis(i).Name)).Kp});
%                 prmlimited{end+1} = tis(i).Name; %#ok<*AGROW>
%                 % TODO correct for permeability limited tissue
%             end
        end
        
        % Add parameters for permeability limited tissues
%         if ~isempty(prmlimited)
%             for i = 1:length(prmlimited)
%                 prm = v_info(strcmpi(v_info(:,1), ...
%                     ['fut_' prmlimited{i} tag]),:);
%                 if ~isempty(prm)
%                     v = addVariantContent(v,{'parameter', ...
%                         ['fut_' prmlimited{i} tag], ...
%                         'Value',kp_prm.Kp(strcmpi({kp_prm.Kp.Tissue}, ...
%                         prmlimited{i})).fut_ic}); 
%                 end
%                 prm = v_info(strcmpi(v_info(:,1), ...
%                     ['Fec_' prmlimited{i}]),:);
%                 if ~isempty(prm) 
%                     v = addVariantContent(v,{'parameter', ...
%                         ['Fec_' prmlimited{i}], ...
%                         'Value',phys.Tissues(strcmpi({phys.Tissues.Name}, ...
%                         prmlimited{i})).Fvec});
%                 end
%             end
%             f_pH7 = smtb.pbpk.ionfrac(inputs.Acidic_pKa,inputs.Acidic_pKa, ...
%                 7);
%             f_pH74 = smtb.pbpk.ionfrac(inputs.Acidic_pKa,inputs.Acidic_pKa, ...
%                 7.4);
%             prm = v_info(strcmpi(v_info(:,1),['PS' tag]),:);     
%             if ~isempty(prm), v = addVariantContent(v,{'parameter', ...
%                     ['PS' tag],         'Value',    0}); end
%             prm = v_info(strcmpi(v_info(:,1),['f_un_pH7' tag]),:); 
%             if ~isempty(prm), v = addVariantContent(v,{'parameter', ...
%                     ['f_un_pH7' tag],   'Value',    f_pH7.neutral}); end
%             prm = v_info(strcmpi(v_info(:,1),['f_un_pH74' tag]),:); 
%             if ~isempty(prm), v = addVariantContent(v,{'parameter',  ...
%                     ['f_un_pH74' tag]   'Value',    f_pH74.neutral}); end
%         end
%         
        % Scale Liver CLint
        lvr_wgt   = phys.Tissues(strcmpi({phys.Tissues.Name}, ...
            'Liver')).Volume * phys.Tissues(strcmpi({phys.Tissues.Name}, ...
            'Liver')).Density;
        if isempty(sp_input.CLint_u_Lhr)
            clh = sp_input.CLint * 60 / 1000 * lvr_wgt ...
                / sp_input.CLint_fuinc;
            sp_input.CLint_u_Lhr = clh;
        else, clh = sp_input.CLint_u_Lhr;
        end
%         v = addVariantContent(v,{'parameter',['CLint_Liver' tag], ...
%             'Value',clh});
        vnew.addcontent({'parameter',['CLint_Liver' tag], ...
            'Value',clh});
        
        % Calculate Plasma Liver Clearance (not used in model)
        lbf = phys.Tissues(strcmp({phys.Tissues.Name}, ...
            'Liver')).TissuePerfusion * 3600 / 1000; % units L/hr
        sp_input.CL_hep_pl = sp_input.B2P * lbf * ...
        (clh/(clh + lbf*(sp_input.B2P/sp_input.adjusted_fup)));
        
        % Calculate Glomeruluar Filtration
        gfr = [phys.Tissues.GlomerularFiltrationRate] * 3600 / 1000; 
        %GFR converted to L/hr
        switch lower(strrep(sp_input.RenalFiltration,' ',''))
            %  case 'fup*gfr', clr = sp_input.fup * gfr / 
            % getVariantVal(v,['fup' tag]);      
            % Need to divide through by the value of fup 
            % used to convert to CLint
            case 'fup*gfr', clr = sp_input.fup * gfr ...
                    / smtb.simbio.getVariantContent(vnew,['fup' tag]);   
                % Need to divide through by the value of fup...
                % used to convert to CLint
                sp_input.CL_ren_pl = gfr * sp_input.fup;                                   
                % Calculate total renal clearance (not used directly 
                % in model)
%             case 'gfr',     clr = gfr / getVariantVal(v,['fup' tag]);                    
% Need to divide through by the value of fup used to convert to CLint
            case 'gfr',     clr = gfr / getVariantContent(v,['fup' tag]);                  
                % Need to divide through by the value of fup used to 
                % convert to CLint
                sp_input.CL_ren_pl = gfr;                                                  
                % Calculate total renal clearance 
                % (not used directly in model)
            case 'off',     clr = 0;
                sp_input.CL_ren_pl = 0;
            case 'fraction*kidneybloodflow'
                warning(['Fraction*Kidney Blood Flow method not' ...
                    ' supported, renal filtration will be set to off'])
                sp_input.RenalFiltration = 'off';
            case 'userdefined'
                clr = sp_input.CL_ren_pl /...
                getVariantContent(v,['fup' tag]); 
                
            otherwise, error('Unrecognized renal filtration method "%s"' ...
                    ,sp_input.RenalFiltration)
        end
%         addcontent(v,{'parameter',['CLint_Kidney' tag],'Value',clr});
%         v = addVariantContent(v,{'parameter',['CLint_Kidney' tag], ...
%             'Value',clr});
        vnew.addcontent({'parameter',['CLint_Kidney' tag], ...
            'Value',clr});
  
        sp_input.CL_tot_pl = ((sp_input.CL_hep_pl + (sp_input.CL_ren_pl)));
        sp_input.CL_tot_bl = sp_input.CL_tot_pl / sp_input.B2P;
        sp_input.CL_units = 'L/hr';
        sp_input.CLh_pctLBF = ((sp_input.CL_hep_pl / sp_input.B2P) ...
            / lbf) * 100;
        
        if isempty(sp_input.CLint)
            sp_input.CLint = (clh * sp_input.CLint_fuinc) / ((60 ...
                / 1000) * lvr_wgt);
        end
        
%     case 'compartmental'
%         v = addVariantContent(v,{'compartment','central','capacity',1});
%         v = addVariantContent(v,{'compartment','peripheral1', ...
%             'capacity',1});
%         v = addVariantContent(v,{'compartment','peripheral2', ...
%             'capacity',1});
%         v = addVariantContent(v,{'parameter','k12','value',0});
%         v = addVariantContent(v,{'parameter','k13','value',0});
%         v = addVariantContent(v,{'parameter','k21','value',0});
%         v = addVariantContent(v,{'parameter','k31','value',0});
%         v = addVariantContent(v,{'parameter','CL','value',0});
%         v = addVariantContent(v,{'parameter','Q_Gut','Value',1}); 
%         % assume G+ makes Q_Gut = 1 in compartmental model
%         v = addVariantContent(v,{'parameter','Q_Liver','Value', ...
%             tis(strcmp({tis.Name},'Liver')).TissuePerfusion});
%         v = addVariantContent(v,{'parameter','Liver','Value', ...
%             tis(strcmp({tis.Name},'Liver')).Volume});
    otherwise, error('Unrecognized distribution model: %s', ...
            inputs.Distribution)
end

end

