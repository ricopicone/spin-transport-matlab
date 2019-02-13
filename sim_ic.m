% --------------------------------------------------------------
function u0 = sim_ic(rr)

global rho_1_langevin rho_2_langevin rho_3_langevin sim rmax

switch sim
    case 'rugar1'
%         u0 = ...
%             [...
%                 rho_1_langevin(rr); ...
%                 rho_2_langevin(rr); ...
%                 rho_3_langevin(rr) ...
%                     * hole(rr,-rmax/18,rmax/16,6e-1) * hole(rr,rmax/18,rmax/16,6e-1) ...
%                     * hole(rr,-2*rmax/16,rmax/16,6e-1) * hole(rr,2*rmax/16,rmax/16,6e-1);
%             ];
        u0 = ...
            [...
                rho_1_langevin(rr); ...
                rho_2_langevin(rr); ...
                rho_3_langevin(rr) ...
                    * -hole_frac(rr,0,rmax/20,1e-1,.01);
            ];
%         u0 = ...
%             [...
%                 rho_1_langevin(rr); ...
%                 rho_2_langevin(rr); ...
%                 rho_3_langevin(rr);
%             ]; 
    case 'performed'
%         w = 4*rmax/20;
%         r0 = -rmax + w/5;
%         slope = 2e1;
%         u0 = ...
%             [...
%                 rho_1_langevin(rr); ...
%                 rho_2_langevin(rr); ...
%                 rho_3_langevin(rr) * ( 1/2 + hole(rr,r0,w,slope)/2 );
%             ];
        u0 = ...
            [...
                rho_1_langevin(rr); ...
                rho_2_langevin(rr); ...
                rho_3_langevin(rr) * ( 1 - exp( - ( rr + rmax )^2/1e0 ) );
            ];
    case 'equilibrium'
%         u0 = ...
%             [...
%                 rho_1_langevin(rr); ...
%                 rho_2_langevin(rr); ...
%                 rho_3_langevin(rr) ...
%             ];
        u0 = ...
            [...
                0; ...
                1 * ( exp( - ( rr )^2/1e-1 ) ); ...
                1 * ( exp( - ( rr )^2/1e-1 ) ) ...
            ];
    otherwise
        u0 = ...
            [...
                rho_1_langevin(rr); ...
                rho_2_langevin(rr); ...
                rho_3_langevin(rr);
            ];
end

end
