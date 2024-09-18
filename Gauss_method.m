%% Power flow analysis using Gauss method
Sbase = 100; Vbase = 230;
%% Lines data: Shunt Y is nominated by Yshunti_j
% Line 1-2:
R1_2 = 0.01008; X1_2 = 0.05040i; Yshunt1_2 = 0.05125i;
Y1_2 = 1/(R1_2 + X1_2);
% Line 1-3:
R1_3 = 0.00744; X1_3 = 0.03720i; Yshunt1_3 = 0.03975i;
Y1_3 = 1/(R1_3 + X1_3);
% Line 2-4:
R2_4 = 0.00744; X2_4 = 0.03720i; Yshunt2_4 = 0.03875i;
Y2_4 = 1/(R2_4 + X2_4);
% Line 3_4:
R3_4 = 0.01272; X3_4 = 0.06360i; Yshunt3_4 = 0.06375i;
Y3_4 = 1/(R3_4 + X3_4);
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Buses data:
% Bus 1: Slack bus
Pload_1 = 50; Qload_1 = 30.99; V_1 = 1;
Ploadpu_1 = Pload_1/Sbase; Qloadpu_1 = Qload_1/Sbase;
% Bus 2: Load bus (inductive)
Pload_2 = 170; Qload_2 = 105.35; V_2 = 1;
Ploadpu_2 = -Pload_2/Sbase; Qloadpu_2 = -Qload_2/Sbase;
% Bus 3: Load bus (inductive)
Pload_3 = 200; Qload_3 = 123.94; V_3 = 1;
Ploadpu_3 = -Pload_3/Sbase; Qloadpu_3 = -Qload_3/Sbase;
% Bus 4: Voltage controlled
Pgen_4 = 318; Pload_4 = 80; Qload_4 = 49.58; V_4 = 1.02;
Ploadpu_4 = (Pgen_4 - Pload_4)/Sbase; Q_4 = 0;
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Admittance bus:
Y11 = Y1_2 + Y1_3 + Yshunt1_2 + Yshunt1_3; Y12 = -Y1_2; Y13 = -Y1_3; Y14 = 0;
Y21 = -Y1_2; Y22 = Y1_2 + Y2_4 + Yshunt1_2 + Yshunt2_4; Y23 = 0; Y24 = -Y2_4;
Y31 = -Y1_3; Y32 = 0; Y33 = Y1_3 + Y3_4 + Yshunt1_3 + Yshunt3_4; Y34 = -Y3_4;
Y41 = 0; Y42 = -Y2_4; Y43 = -Y3_4; Y44 = Y2_4 + Y3_4 + Yshunt2_4 + Yshunt3_4;
Y = [ Y11 Y12 Y13 Y14 ; Y21 Y22 Y23 Y24 ; Y31 Y32 Y33 Y34 ; Y41 Y42 Y43 Y44];
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Define Vitertion array which contains the current iteration of each bus
%% voltage, which starts with the initial conditions as the following:
Viteration = [ V_1 ; V_2 ; V_3 ; V_4 ];
%% Define number of iterations: (15 iterations was selected to be the
%% required number of iterations till steady state results by trial & error)
Niteration = 15;
%% Bus voltage calculations based on the selected number of iterations:
for i = 1:Niteration
    for n = 2:4
         if n == 2
             Viteration(2) = (1/Y(2,2))*(((Ploadpu_2 - Qloadpu_2*1i)/(conj(Viteration(2)))) - Y(2,1)*Viteration(1) - Y(2,3)*Viteration(3) - Y(2,4)*Viteration(4));
         elseif n == 3
             Viteration(3) = (1/Y(3,3))*(((Ploadpu_3 - Qloadpu_3*1i)/(conj(Viteration(3)))) - Y(3,1)*Viteration(1) - Y(3,2)*Viteration(2) - Y(3,4)*Viteration(4));
         elseif n == 4
             Q_4 = -imag((conj(Viteration(4))*(Y(4,1)*Viteration(1) + Y(4,2)*Viteration(2) + Y(4,3)*Viteration(3) + Y(4,4)*Viteration(4))));
             Viteration(4) = (1/Y(4,4))*(((Ploadpu_4 - Q_4*1i)/(conj(Viteration(4)))) - Y(4,1)*Viteration(1) - Y(4,2)*Viteration(2) - Y(4,3)*Viteration(3));
             Viteration(4) = V_4*(cos(angle(Viteration(4))) + sin(angle(Viteration(4)))*1i);
        end
    end
end 
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Slack bus power calculation:
Pslack_inj = real(V_1*conj((Y(1,1)*V_1 + Y(1,2)*Viteration(2) + Y(1,3)*Viteration(3) + Y(1,4)*Viteration(4))));
Qslack_inj = imag(V_1*conj((Y(1,1)*V_1 + Y(1,2)*Viteration(2) + Y(1,3)*Viteration(3) + Y(1,4)*Viteration(4))));
Pslack_gen = Pslack_inj + Ploadpu_1;
Qslack_gen = Qslack_inj + Qloadpu_1;
Sslack_inj = Pslack_inj + Qslack_inj*1i;
Sslack_gen = Pslack_gen + Qslack_gen*1i;
%% Line flow & power losses:
% Line 1-2:
S1_2 = V_1*(conj((V_1 - Viteration(2))*Y1_2 + V_1*Yshunt1_2));
S2_1 = Viteration(2)*(conj((Viteration(2) - V_1)*Y1_2 + Viteration(2)*Yshunt1_2));
Slosses1_2 = S1_2 + S2_1;
% Line 1-3:
S1_3 = V_1*(conj((V_1 - Viteration(3))*Y1_3 + V_1*Yshunt1_3));
S3_1 = Viteration(3)*(conj((Viteration(3) - V_1)*Y1_3 + Viteration(3)*Yshunt1_3));
Slosses1_3 = S1_3 + S3_1;
% Line 2-4:
S2_4 = Viteration(2)*(conj((Viteration(2) - Viteration(4))*Y2_4 + Viteration(2)*Yshunt2_4));
S4_2 = Viteration(4)*(conj((Viteration(4) - Viteration(2))*Y2_4 + Viteration(4)*Yshunt2_4));
Slosses2_4 = S2_4 + S4_2;
% Line 3-4:
S3_4 = Viteration(3)*(conj((Viteration(3) - Viteration(4))*Y3_4 + Viteration(3)*Yshunt3_4));
S4_3 = Viteration(4)*(conj((Viteration(4) - Viteration(3))*Y3_4 + Viteration(4)*Yshunt3_4));
Slosses3_4 = S3_4 + S4_3;
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Results display:

%Bus voltages display:
Bus_number = ["Bus 1";"Bus 2"; "Bus 3"; "Bus 4"];
Bus_voltage_in_kV =Vbase*Viteration;
table(Bus_number,Bus_voltage_in_kV)

%Slack bus power display:
Slack_bus_data = ["Injected power in MVA";"Generated power in MVA"];
Slack_bus_power = Sbase*[Sslack_inj; Sslack_gen];
table(Slack_bus,Slack_bus_power)

%Line flows and lines losses display:
Line =["Line 1-2";"Line 1-3";"Line 2-4";"Line 3-4";];
Line_flow = Sbase*[S1_2; S1_3; S4_2; S4_3];
Direction = ["Bus 1 to Bus 2";"Bus 1 to Bus 3";"Bus 4 to Bus 3";"Bus 4 to Bus 2";];
Line_losses = Sbase*[Slosses1_2; Slosses1_3; Slosses2_4; Slosses3_4];
table(Line, Line_flow,Direction ,Line_losses)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
