function counter()
% Simulate a two-state model of gene expression
import Gillespie.*

%% Reaction network:
%%% DNA_0
%   1. TF activation:      DNA_0 + Ara              --k_BAD_on-->                   DNA_0_BAD
%   2. TF deactivation:    DNA_0_BAD                --k_BAD_off-->                  DNA_0 + Ara

%   3. transcription:      DNA_0_BAD                --kM_BAD_flp-->                 DNA_0_BAD + mRNA_flp

%   4. translation:        mRNA_flp                 --kP_Flp-->                     mRNA_flp + Flp

%   ((5. recombination:      DNA_0 + Flp                --kR_Flp-->                     DNA_1 + Flp)) ignored for simplicity

%   6. recombination:      DNA_0_BAD + Flp          --kR_Flp-->                     DNA_1_BAD + Flp

%   7. mRNA_flp decay:     mRNA_flp                 --gM_flp-->                     0
%   8. Flp decay:          Flp                      --gP_Flp-->                     0


%%% DNA_1
%   9. TF activation:      DNA_1 + Ara              --k_BAD_on-->                   DNA_1_BAD  (same rate as for DNA_0)
%   10. TF deactivation:   DNA_1_BAD                --k_BAD_off-->                  DNA_1 + Ara

%   11. transcription:     DNA_1_BAD                --kM_BAD_gfp_iptg-->            DNA_1_BAD + mRNA_gfp + mRNA_iptg

%   12. translation:       mRNA_gfp                 --kP_GFP-->                     mRNA_gfp + GFP
%   13. translation:       mRNA_iptg                --kP_IPTG-->                    mRNA_iptg + IPTG

%   ((14. TF activation:      DNA_1 + IPTG              --k_A1lacO_on-->                DNA_1_A1lacO))
%   ((15. TF deactivation:    DNA_1_A1lacO              --k_A1lacO_off-->               DNA_1 + IPTG))

%   16. TF activation:     DNA_1_BAD + IPTG         --k_A1lacO_on-->                DNA_1_BAD_A1lacO
%   17. TF deactivation:   DNA_1_BAD_A1lacO         --k_A1lacO_off-->               DNA_1_BAD + IPTG

%   ((18. TF activation:     DNA_1_A1lacO + Aba         --k_BAD_on-->                   DNA_1_BAD_A1lacO))
%   ((19. TF deactivation:   DNA_1_BAD_A1lacO           --k_BAD_off-->                  DNA_1_A1lacO + Aba))

%   ((20. transcription:     DNA_1_A1lacO               --kM_IPTG_flp-->                DNA_1_A1lacO + mRNA_flp))

%   21. transcription:     DNA_1_BAD_A1lacO         --kM_BAD_IPTG_gfp_iptg_flp-->   DNA_1_BAD_A1lacO + mRNA_gfp + mRNA_iptg + mRNA_flp

%   ((22. recombination:     DNA_1 + Flp                --kR_Flp-->                     DNA_0 + Flp))
%   ((23. recombination:     DNA_1_BAD + Flp            --kR_Flp-->                     DNA_0_BAD + Flp))
%   ((24. recombination:     DNA_1_A1lacO + Flp         --kR_Flp-->                     DNA_0_A1lacO + Flp))
%   25. recombination:     DNA_1_BAD_A1lacO + Flp   --kR_Flp-->                     DNA_0_BAD_A1lacO + Flp

%   26. mRNA_gfp decay:    mRNA_gfp                 --gM_gfp-->                     0
%   27. mRNA_iptg decay:   mRNA_iptg                --gM_iptg-->                    0
%   28. GFP decay:         GFP                      --gP_GFP-->                     0
%   29. IPTG decay:        IPTG                     --gP_IPTG-->                    0
%   30. Ara decay:         Ara                      --gP_Ara-->                     0







%% Rate constants
p.k_BAD_on = 0.0002; % sec^-1 -> calculated by hand
p.k_BAD_off = 0.01; %sec^-1 -> given by kobi
p.kM_BAD_flp = 0.1; 
p.kP_Flp = 0.1;
p.kR_Flp = 0.1;
p.gM_flp =0.1;
p.gP_Flp  = 0.1;
p.kM_BAD_gfp_iptg = 0.1;
p.kP_GFP = 0.1;
p.kP_IPTG = 0.1;
p.k_A1lacO_on = 0.1;
p.k_A1lacO_off = 0.1;
p.kM_BAD_IPTG_gfp_iptg_flp = 0.1;
p.gM_gfp = 0.1;
p.gM_iptg = 0.1;
p.gP_GFP = 0.1;
p.gP_IPTG = 0.1;
p.gP_Ara = 0.1201;
%% Initial state
tspan = [0, 10000]; %seconds
x0    = [100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];     %DNA_0, Ara, DNA_O_BAD, mRNA_flp, Flp, DNA_1, DNA_1_BAD, mRNA_gfp, mRNA_iptg, GFP, IPTG, DNA_1_BAD_A1lacO, DNA_0_BAD_A1lacO

%% Specify reaction network
pfun = @propensities_2state;
                %DNA_0, Ara, DNA_O_BAD, mRNA_flp, Flp, DNA_1, DNA_1_BAD, mRNA_gfp, mRNA_iptg, GFP, IPTG, DNA_1_BAD_A1lacO, DNA_0_BAD_A1lacO
stoich_matrix = [-1     -1      1          0       0     0        0         0          0       0     0          0               0            %DNA_0 + Ara              --k_BAD_on-->                   DNA_0_BAD
                  1      1     -1          0       0     0        0         0          0       0     0          0               0            %DNA_0_BAD                --k_BAD_off-->                  DNA_0 + Ara
                  0      0      0          1       0     0        0         0          0       0     0          0               0            %DNA_0_BAD                --kM_BAD_flp-->                 DNA_0_BAD + mRNA_flp
                  0      0      0          0       1     0        0         0          0       0     0          0               0            %mRNA_flp                 --kP_Flp-->                     mRNA_flp + Flp
                  0      0     -1          0       0     0        1         0          0       0     0          0               0            %DNA_0_BAD + Flp          --kR_Flp-->                     DNA_1_BAD + Flp
                  0      0      0         -1       0     0        0         0          0       0     0          0               0            %mRNA_flp                 --gM_flp-->                     0
                  0      0      0          0      -1     0        0         0          0       0     0          0               0            %Flp                      --gP_Flp-->                     0
                  0     -1      0          0       0    -1        1         0          0       0     0          0               0            %DNA_1 + Ara              --k_BAD_on-->                   DNA_1_BAD  (same rate as for DNA_0)
                  0      1      0          0       0     1       -1         0          0       0     0          0               0            %DNA_1_BAD                --k_BAD_off-->                  DNA_1 + Ara
                  0      0      0          0       0     0        0         1          1       0     0          0               0            %DNA_1_BAD                --kM_BAD_gfp_iptg-->            DNA_1_BAD + mRNA_gfp + mRNA_iptg
                  0      0      0          0       0     0        0         0          0       1     0          0               0            %mRNA_gfp                 --kP_GFP-->                     mRNA_gfp + GFP
                  0      0      0          0       0     0        0         0          0       0     1          0               0            %mRNA_iptg                --kP_IPTG-->                    mRNA_iptg + IPTG
                  0      0      0          0       0     0       -1         0          0       0    -1          1               0            %DNA_1_BAD + IPTG         --k_A1lacO_on-->                DNA_1_BAD_A1lacO
                  0      0      0          0       0     0        1         0          0       0     1         -1               0            %DNA_1_BAD_A1lacO         --k_A1lacO_off-->               DNA_1_BAD + IPTG
                  0      0      0          1       0     0        0         1          1       0     0          0               0            %DNA_1_BAD_A1lacO         --kM_BAD_IPTG_gfp_iptg_flp-->   DNA_1_BAD_A1lacO + mRNA_gfp + mRNA_iptg + mRNA_flp
                  0      0      0          0       0     0        0         0          0       0     0         -1               1            %DNA_1_BAD_A1lacO + Flp   --kR_Flp-->                     DNA_0_BAD_A1lacO + Flp
                  0      0      0          0       0     0        0        -1          0       0     0          0               0            %mRNA_gfp                 --gM_gfp-->                     0
                  0      0      0          0       0     0        0         0         -1       0     0          0               0            %mRNA_iptg                --gM_iptg-->                    0
                  0      0      0          0       0     0        0         0          0      -1     0          0               0            %GFP                      --gP_GFP-->                     0
                  0      0      0          0       0     0        0         0          0       0    -1          0               0            %IPTG                     --gP_IPTG-->                    0
                  0     -1      0          0       0     0        0         0          0       0     0          0               0];          %Ara                      --gP_Ara-->                     0
                

%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);

%% Plot time course
figure();
stairs(t,x); 
set(gca,'XLim',tspan);
xlabel('time (s)');
ylabel('molecules');
legend({'DNA_0', 'Ara', 'DNA_O_BAD', 'mRNA_flp', 'Flp', 'DNA_1', 'DNA_1_BAD', 'mRNA_gfp', 'mRNA_iptg', 'GFP', 'IPTG', 'DNA_1_BAD_A1lacO', 'DNA_0_BAD_A1lacO'});

end


function a = propensities_2state(x, p)
% Return reaction propensities given current state x
DNA_0    = x(1);
Ara = x(2);
mRNA_flp = x(4);

a = [p.k_BAD_on*DNA_0;            %transcription
     p.k_BAD_on*Ara;       %translation
     p.k_BAD_on*DNA_0;       %mRNA decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp;   %protein decay
     p.k_BAD_on*mRNA_flp];   %protein decay 
end
