clear;clc;
disp('Note:')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%An open-source code to generate carbon nanotube/graphene junctions

%createCNT2D.m is the main function to generate CNT/graphene structure, which only allows carbon nanotbue to rotate along Z direction.
%fun.m and JoiningCNT.m are subfunctions for createCNT2D.m to generate least squares euqations and solve the euqations.

%createCNT3D.m is the main function to generate CNT/graphene structure, which allows carbon nanotbue to rotate along 3 directions.
%fun2.m and JoiningCNT2.m are subfunctions for createCNT3D.m to generate least squares euqations and solve the euqations.

%12/29/2017

%% In creatCNTG_V1.01, sandwich structure is allowed to generate by setting sandwich=1.
%02/28/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('An open-source code to generate carbon nanotube/graphene junctions')
disp(' ')
disp('createCNT2D.m is the main function to generate CNT/graphene structure, which only allows carbon nanotbue to rotate along Z direction.')
disp('fun.m and JoiningCNT.m are subfunctions for createCNT2D.m to generate least squares euqations and solve the euqations.')
disp(' ')
disp('createCNT3D.m is the main function to generate CNT/graphene structure, which allows carbon nanotbue to rotate along 3 directions.')
disp('fun2.m and JoiningCNT2.m are subfunctions for createCNT3D.m to generate least squares euqations and solve the euqations.')
disp(' ')
disp('12/29/2017')

disp('In creatCNTG_V1.01, sandwich structure is allowed to generate by setting sandwich=1.')
disp(' ')
disp('12/29/2017')
