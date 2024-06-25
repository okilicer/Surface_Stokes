function [point_val, hess_ak_elem1, hess_ak_elem2, hess_ak_elem3, grad_ak_elem] = get_grad_ak(point, elem_number, surf_deg, node, elem, surfacedata)

pt2 = point(2);
pt3 = point(3); 


compx1 = 0;    %s Derivative of x-component of poly interp of geometry
compy1 = 0;    %s Derivative of y-component of poly interp of geometry
compz1 = 0;    %s Derivative of z-component of poly interp of geometry

compx2 = 0;    %t Derivative of x-component of poly interp of geometry
compy2 = 0;    %t Derivative of y-component of poly interp of geometry
compz2 = 0;    %t Derivative of z-component of poly interp of geometry 
i = elem_number;

point_comp1 = 0;
point_comp2 = 0;
point_comp3 = 0;

hess_compx1 = 0;
hess_compy1 = 0;
hess_compz1 = 0;
hess_compx2 = 0;
hess_compy2 = 0;
hess_compz2 = 0;
hess_compx3 = 0;
hess_compy3 = 0;
hess_compz3 = 0;
hess_compx4 = 0;
hess_compy4 = 0;
hess_compz4 = 0;

switch surf_deg
    case 1 
        matrix = [-1 -1; 
                   1 0; 
                   0 1];   
        hessian_matrix = [0 0 0 0; 0 0 0 0; 0 0 0 0];
        point_matrix = [1-pt2-pt3;
                        pt2;
                        pt3];
   
        % 6 1 2 order 1-x-y, x, y 
        % %Loop over interpolation points of element
        for m=1:3
            %Derivatives of interpolation of geometry for building metric
            compx1 = compx1 + matrix(m, 1) * node(elem(i, m), 1);
            compy1 = compy1 + matrix(m, 1) * node(elem(i, m), 2);
            compz1 = compz1 + matrix(m, 1) * node(elem(i, m), 3);
            
            compx2 = compx2 + matrix(m, 2) * node(elem(i, m), 1);
            compy2 = compy2 + matrix(m, 2) * node(elem(i, m), 2);            
            compz2 = compz2 + matrix(m, 2) * node(elem(i, m), 3); 

            point_comp1 = point_comp1 + point_matrix(m, 1) * node(elem(i, m), 1);
            point_comp2 = point_comp2 + point_matrix(m, 1) * node(elem(i, m), 2);
            point_comp3 = point_comp3 + point_matrix(m, 1) * node(elem(i, m), 3);
            

        end

        grad_ak_elem = [compx1, compx2; compy1 compy2; compz1 compz2]; 
        hess_ak_elem1 = [0, 0; 0, 0];
        hess_ak_elem2 = [0, 0; 0, 0];
        hess_ak_elem3 = [0, 0; 0, 0];
        point_val = [point_comp1, point_comp2, point_comp3];

    case 2
        % 6 1 2 new points between 6-1, 6-2 and 1-2
        %the functions are
        % 1-3*x-3*y+4*x*y+2*x^2+2*y^2, 2*x^2-x, 2*y^2-y,
        % 4*x-4*x^2-4*x*y, 4*y-4*x*y-4*y^2, 4*x*y
        % the new points are 1/2*node(elem(1, 1), :) + 1/2*node(elem(1, 2), :)
        % 1/2*node(elem(1, 1), :) + 1/2*node(elem(1, 3), :)
        %and 1/2*node(elem(1, 2), :) + 1/2*node(elem(1, 3), :)
        point_matrix = [1-3*pt2-3*pt3+4*pt2*pt3+2*pt2^2+2*pt3^2;
                        2*pt2^2-pt2;
                        2*pt3^2-pt3;
                        4*pt2-4*pt2^2-4*pt2*pt3;
                        4*pt3-4*pt2*pt3-4*pt3^2;
                         4*pt2*pt3];
        
        
        matrix = [-3+4*pt3+4*pt2 -3+4*pt2+4*pt3; 
                   4*pt2-1 0; 0 4*pt3-1;                  
                   4-8*pt2-4*pt3 -4*pt2;
                   -4*pt3 4-8*pt3-4*pt2;
                   4*pt3 4*pt2];
        hess_matrix = [4 4 4 4; 4 0 0 0; 0 0 0 4;-8 -4 -4 0; 0 -4 -4 -8; 0 4 4 0];
        % 6 1 2 order 1-x-y, x, y 
        % %Loop over interpolation points of element
        node_temp = [node(elem(i, 1), :); node(elem(i, 2), :); node(elem(i, 3), :)];
        node1 = 1/2*node(elem(i, 1), :) + 1/2*node(elem(i, 2), :);
%         node1 = 1.1*node1 /norm(node1, 2);
        node1 = surfacedata.project(node1); 
        node2 = 1/2*node(elem(i, 1), :) + 1/2*node(elem(i, 3), :);
%         node2 = 1.1*node2 /norm(node2, 2);
        node2 = surfacedata.project(node2); 
        node3 = 1/2*node(elem(i, 2), :) + 1/2*node(elem(i, 3), :);
%         node3 = 1.1*node3 /norm(node3, 2);
        node3 = surfacedata.project(node3); 
        node_temp(4:6, :) = [node1;
                             node2;
                             node3];


        for m=1:6
            %Derivatives of interpolation of geometry for building metric
            compx1 = compx1 + matrix(m, 1) * node_temp(m, 1);
            compy1 = compy1 + matrix(m, 1) * node_temp(m, 2);
            compz1 = compz1 + matrix(m, 1) * node_temp(m, 3);
            
            compx2 = compx2 + matrix(m, 2) * node_temp(m, 1);
            compy2 = compy2 + matrix(m, 2) * node_temp(m, 2);            
            compz2 = compz2 + matrix(m, 2) * node_temp(m, 3);

            hess_compx1 = hess_compx1 + hess_matrix(m, 1) * node_temp(m, 1);
            hess_compy1 = hess_compy1 + hess_matrix(m, 1) * node_temp(m, 2);
            hess_compz1 = hess_compz1 + hess_matrix(m, 1) * node_temp(m, 3);

            hess_compx2 = hess_compx2 + hess_matrix(m, 2) * node_temp(m, 1);
            hess_compy2 = hess_compy2 + hess_matrix(m, 2) * node_temp(m, 2);
            hess_compz2 = hess_compz2 + hess_matrix(m, 2) * node_temp(m, 3);

            hess_compx3 = hess_compx3 + hess_matrix(m, 3) * node_temp(m, 1);
            hess_compy3 = hess_compy3 + hess_matrix(m, 3) * node_temp(m, 2);
            hess_compz3 = hess_compz3 + hess_matrix(m, 3) * node_temp(m, 3);

            hess_compx4 = hess_compx4 + hess_matrix(m, 4) * node_temp(m, 1);
            hess_compy4 = hess_compy4 + hess_matrix(m, 4) * node_temp(m, 2);
            hess_compz4 = hess_compz4 + hess_matrix(m, 4) * node_temp(m, 3);

            point_comp1 = point_comp1 + point_matrix(m, 1) * node_temp(m, 1);
            point_comp2 = point_comp2 + point_matrix(m, 1) * node_temp(m, 2);
            point_comp3 = point_comp3 + point_matrix(m, 1) * node_temp(m, 3);
        end

        grad_ak_elem = [compx1, compx2; compy1 compy2; compz1 compz2]; 
        hess_ak_elem1 = [hess_compx1, hess_compx2; hess_compx3,hess_compx4];
        hess_ak_elem2 = [hess_compy1, hess_compy2; hess_compy3,hess_compy4];
        hess_ak_elem3 = [hess_compz1, hess_compz2; hess_compz3,hess_compz4];
        point_val = [point_comp1, point_comp2, point_comp3];



    case 3
        point_matrix = [-9*pt2^3/2 - 27*pt2^2*pt3/2 + 9*pt2^2 - 27*pt2*pt3^2/2 + 18*pt2*pt3 - 11*pt2/2 - 9*pt3^3/2 + 9*pt3^2 - 11*pt3/2 + 1;
                        9*pt2^3/2 - 9*pt2^2/2 + pt2;
                        9*pt3^3/2 - 9*pt3^2/2 + pt3;
                        27*pt2^3/2 + 27*pt2^2*pt3 - 45*pt2^2/2 + 27*pt2*pt3^2/2 - 45*pt2*pt3/2 + 9*pt2;
                        -27*pt2^3/2 - 27*pt2^2*pt3/2 + 18*pt2^2 + 9*pt2*pt3/2 - 9*pt2/2;
                         27*pt2^2*pt3/2 - 9*pt2*pt3/2;
                        27*pt2*pt3^2/2 - 9*pt2*pt3/2;
                        27*pt2^2*pt3/2 + 27*pt2*pt3^2 - 45*pt2*pt3/2 + 27*pt3^3/2 - 45*pt3^2/2 + 9*pt3;
                        -27*pt2*pt3^2/2 + 9*pt2*pt3/2 - 27*pt3^3/2 + 18*pt3^2 - 9*pt3/2;
                        -27*pt2^2*pt3 - 27*pt2*pt3^2 + 27*pt2*pt3];
        matrix = [ 18*pt2 + 18*pt3 - 27*pt2*pt3 - (27*pt2^2)/2 - (27*pt3^2)/2 - 11/2,  18*pt2 + 18*pt3 - 27*pt2*pt3 - (27*pt2^2)/2 - (27*pt3^2)/2 - 11/2;
                                  (27*pt2^2)/2 - 9*pt2 + 1,                                                      0;
                                                     0,                                   (27*pt3^2)/2 - 9*pt3 + 1;
                            (81*pt2^2)/2 + 54*pt2*pt3 - 45*pt2 + (27*pt3^2)/2 - (45*pt3)/2 + 9,                             27*pt2*pt3 - (45*pt2)/2 + 27*pt2^2;
                            36*pt2 + (9*pt3)/2 - 27*pt2*pt3 - (81*pt2^2)/2 - 9/2,                                   (9*pt2)/2 - (27*pt2^2)/2;
                                      27*pt2*pt3 - (9*pt3)/2,                                   (27*pt2^2)/2 - (9*pt2)/2;
                                  (27*pt3^2)/2 - (9*pt3)/2,                                       27*pt2*pt3 - (9*pt2)/2;
                             27*pt2*pt3 - (45*pt3)/2 + 27*pt3^2, (27*pt2^2)/2 + 54*pt2*pt3 - (45*pt2)/2 + (81*pt3^2)/2 - 45*pt3 + 9;
                                  (9*pt3)/2 - (27*pt3^2)/2,             (9*pt2)/2 + 36*pt3 - 27*pt2*pt3 - (81*pt3^2)/2 - 9/2;
                                27*pt3 - 54*pt2*pt3 - 27*pt3^2,                                 27*pt2 - 54*pt2*pt3 - 27*pt2^2];
        hess_matrix = [18 - 27*pt3 - 27*pt2,   18 - 27*pt3 - 27*pt2,   18 - 27*pt3 - 27*pt2, 18 - 27*pt3 - 27*pt2;
                    27*pt2 - 9,                  0,                  0,                0;
                    0,                  0,                  0,         27*pt3 - 9;
                    81*pt2 + 54*pt3 - 45, 54*pt2 + 27*pt3 - 45/2, 54*pt2 + 27*pt3 - 45/2,             27*pt2;
                    36 - 27*pt3 - 81*pt2,         9/2 - 27*pt2,         9/2 - 27*pt2,                0;
                    27*pt3,         27*pt2 - 9/2,         27*pt2 - 9/2,                0;
                     0,         27*pt3 - 9/2,         27*pt3 - 9/2,             27*pt2;
                    27*pt3, 27*pt2 + 54*pt3 - 45/2, 27*pt2 + 54*pt3 - 45/2, 54*pt2 + 81*pt3 - 45;
                    0,         9/2 - 27*pt3,         9/2 - 27*pt3, 36 - 81*pt3 - 27*pt2;
                    -54*pt3,   27 - 54*pt3 - 54*pt2,   27 - 54*pt3 - 54*pt2,            -54*pt2];

         node_temp = [node(elem(i, 1), :); node(elem(i, 2), :); node(elem(i, 3), :)];
         node1 = 2/3*node(elem(i, 1), :) + 1/3*node(elem(i, 2), :);
         %node1 = node1 /norm(node1, 2);
         node1 = surfacedata.project(node1);
         node2 = 1/3*node(elem(i, 1), :) + 2/3*node(elem(i, 2), :);
%          node2 = node2 /norm(node2, 2);
         node2 = surfacedata.project(node2);
         node3 = 2/3*node(elem(i, 2), :) + 1/3*node(elem(i, 3), :);
%          node3 = node3 /norm(node3, 2);
         node3 = surfacedata.project(node3);
         node4 = 1/3*node(elem(i, 2), :) + 2/3*node(elem(i, 3), :);
         node4 = surfacedata.project(node4);
         node5 = 2/3*node(elem(i, 1), :) + 1/3*node(elem(i, 3), :);
         node5 = surfacedata.project(node5);
         node6 = 1/3*node(elem(i, 1), :) + 2/3*node(elem(i, 3), :);
         node6 = surfacedata.project(node6);
         node7 = 1/3*node(elem(i, 1), :) + 1/3*node(elem(i, 2), :) + 1/3*node(elem(i, 3), :);
         node7 = surfacedata.project(node7);
         node_temp(4:10, :) = [node1;
                             node2;
                             node3;
                             node4;
                             node5;
                             node6;
                             node7];
         for m=1:10
            %Derivatives of interpolation of geometry for building metric
            compx1 = compx1 + matrix(m, 1) * node_temp(m, 1);
            compy1 = compy1 + matrix(m, 1) * node_temp(m, 2);
            compz1 = compz1 + matrix(m, 1) * node_temp(m, 3);
            
            compx2 = compx2 + matrix(m, 2) * node_temp(m, 1);
            compy2 = compy2 + matrix(m, 2) * node_temp(m, 2);            
            compz2 = compz2 + matrix(m, 2) * node_temp(m, 3);

            hess_compx1 = hess_compx1 + hess_matrix(m, 1) * node_temp(m, 1);
            hess_compy1 = hess_compy1 + hess_matrix(m, 1) * node_temp(m, 2);
            hess_compz1 = hess_compz1 + hess_matrix(m, 1) * node_temp(m, 3);

            hess_compx2 = hess_compx2 + hess_matrix(m, 2) * node_temp(m, 1);
            hess_compy2 = hess_compy2 + hess_matrix(m, 2) * node_temp(m, 2);
            hess_compz2 = hess_compz2 + hess_matrix(m, 2) * node_temp(m, 3);

            hess_compx3 = hess_compx3 + hess_matrix(m, 3) * node_temp(m, 1);
            hess_compy3 = hess_compy3 + hess_matrix(m, 3) * node_temp(m, 2);
            hess_compz3 = hess_compz3 + hess_matrix(m, 3) * node_temp(m, 3);

            hess_compx4 = hess_compx4 + hess_matrix(m, 4) * node_temp(m, 1);
            hess_compy4 = hess_compy4 + hess_matrix(m, 4) * node_temp(m, 2);
            hess_compz4 = hess_compz4 + hess_matrix(m, 4) * node_temp(m, 3);

            point_comp1 = point_comp1 + point_matrix(m, 1) * node_temp(m, 1);
            point_comp2 = point_comp2 + point_matrix(m, 1) * node_temp(m, 2);
            point_comp3 = point_comp3 + point_matrix(m, 1) * node_temp(m, 3);
        end

        grad_ak_elem = [compx1, compx2; compy1 compy2; compz1 compz2]; 
        hess_ak_elem1 = [hess_compx1, hess_compx2; hess_compx3,hess_compx4];
        hess_ak_elem2 = [hess_compy1, hess_compy2; hess_compy3,hess_compy4];
        hess_ak_elem3 = [hess_compz1, hess_compz2; hess_compz3,hess_compz4];
        point_val = [point_comp1, point_comp2, point_comp3];
    case 4
            point_matrix = [32*pt2^4/3 + 128*pt2^3*pt3/3 - 80*pt2^3/3 + 64*pt2^2*pt3^2 - 80*pt2^2*pt3 + 70*pt2^2/3 + 128*pt2*pt3^3/3 - 80*pt2*pt3^2 + 140*pt2*pt3/3 - 25*pt2/3 + 32*pt3^4/3 - 80*pt3^3/3 + 70*pt3^2/3 - 25*pt3/3 + 1;
                            32*pt2^4/3 - 16*pt2^3 + 22*pt2^2/3 - pt2;
                            32*pt3^4/3 - 16*pt3^3 + 22*pt3^2/3 - pt3;
                            128*pt2^3*pt3/3 - 32*pt2^2*pt3 + 16*pt2*pt3/3;
                            64*pt2^2*pt3^2 - 16*pt2^2*pt3 - 16*pt2*pt3^2 + 4*pt2*pt3;
                            128*pt2*pt3^3/3 - 32*pt2*pt3^2 + 16*pt2*pt3/3;
                            -128*pt2^3*pt3/3 - 128*pt2^2*pt3^2 + 96*pt2^2*pt3 - 128*pt2*pt3^3 + 192*pt2*pt3^2 - 208*pt2*pt3/3 - 128*pt3^4/3 + 96*pt3^3 - 208*pt3^2/3 + 16*pt3;
                            64*pt2^2*pt3^2 - 16*pt2^2*pt3 + 128*pt2*pt3^3 - 144*pt2*pt3^2 + 28*pt2*pt3 + 64*pt3^4 - 128*pt3^3 + 76*pt3^2 - 12*pt3;
                            -128*pt2*pt3^3/3 + 32*pt2*pt3^2 - 16*pt2*pt3/3 - 128*pt3^4/3 + 224*pt3^3/3 - 112*pt3^2/3 + 16*pt3/3;
                            -128*pt2^4/3 - 128*pt2^3*pt3 + 96*pt2^3 - 128*pt2^2*pt3^2 + 192*pt2^2*pt3 - 208*pt2^2/3 - 128*pt2*pt3^3/3 + 96*pt2*pt3^2 - 208*pt2*pt3/3 + 16*pt2;
                            64*pt2^4 + 128*pt2^3*pt3 - 128*pt2^3 + 64*pt2^2*pt3^2 - 144*pt2^2*pt3 + 76*pt2^2 - 16*pt2*pt3^2 + 28*pt2*pt3 - 12*pt2;
                            -128*pt2^4/3 - 128*pt2^3*pt3/3 + 224*pt2^3/3 + 32*pt2^2*pt3 - 112*pt2^2/3 - 16*pt2*pt3/3 + 16*pt2/3;
                            128*pt2^3*pt3 + 256*pt2^2*pt3^2 - 224*pt2^2*pt3 + 128*pt2*pt3^3 - 224*pt2*pt3^2 + 96*pt2*pt3;
                            -128*pt2^3*pt3 - 128*pt2^2*pt3^2 + 160*pt2^2*pt3 + 32*pt2*pt3^2 - 32*pt2*pt3;
                            -128*pt2^2*pt3^2 + 32*pt2^2*pt3 - 128*pt2*pt3^3 + 160*pt2*pt3^2 - 32*pt2*pt3];

            matrix = [ (128*pt2^3)/3 + 128*pt2^2*pt3 - 80*pt2^2 + 128*pt2*pt3^2 - 160*pt2*pt3 + (140*pt2)/3 + (128*pt3^3)/3 - 80*pt3^2 + (140*pt3)/3 - 25/3,  (128*pt2^3)/3 + 128*pt2^2*pt3 - 80*pt2^2 + 128*pt2*pt3^2 - 160*pt2*pt3 + (140*pt2)/3 + (128*pt3^3)/3 - 80*pt3^2 + (140*pt3)/3 - 25/3;
                       (128*pt2^3)/3 - 48*pt2^2 + (44*pt2)/3 - 1,                                                                                                             0;
                       0,                                                                           (128*pt3^3)/3 - 48*pt3^2 + (44*pt3)/3 - 1;
                       128*pt3*pt2^2 - 64*pt3*pt2 + (16*pt3)/3,                                                                               (128*pt2^3)/3 - 32*pt2^2 + (16*pt2)/3;
                       4*pt3 - 32*pt2*pt3 + 128*pt2*pt3^2 - 16*pt3^2,                                                                             4*pt2 - 32*pt2*pt3 + 128*pt2^2*pt3 - 16*pt2^2;
                       (128*pt3^3)/3 - 32*pt3^2 + (16*pt3)/3,                                                                                 128*pt2*pt3^2 - 64*pt2*pt3 + (16*pt2)/3;
                       - 128*pt2^2*pt3 - 256*pt2*pt3^2 + 192*pt2*pt3 - 128*pt3^3 + 192*pt3^2 - (208*pt3)/3, - (128*pt2^3)/3 - 256*pt2^2*pt3 + 96*pt2^2 - 384*pt2*pt3^2 + 384*pt2*pt3 - (208*pt2)/3 - (512*pt3^3)/3 + 288*pt3^2 - (416*pt3)/3 + 16;
                       28*pt3 - 32*pt2*pt3 + 128*pt2*pt3^2 - 144*pt3^2 + 128*pt3^3,                              128*pt2^2*pt3 - 16*pt2^2 + 384*pt2*pt3^2 - 288*pt2*pt3 + 28*pt2 + 256*pt3^3 - 384*pt3^2 + 152*pt3 - 12;
                       - (128*pt3^3)/3 + 32*pt3^2 - (16*pt3)/3,                                      64*pt2*pt3 - (224*pt3)/3 - (16*pt2)/3 - 128*pt2*pt3^2 + 224*pt3^2 - (512*pt3^3)/3 + 16/3;
                    - (512*pt2^3)/3 - 384*pt2^2*pt3 + 288*pt2^2 - 256*pt2*pt3^2 + 384*pt2*pt3 - (416*pt2)/3 - (128*pt3^3)/3 + 96*pt3^2 - (208*pt3)/3 + 16,                                             - 128*pt2^3 - 256*pt2^2*pt3 + 192*pt2^2 - 128*pt2*pt3^2 + 192*pt2*pt3 - (208*pt2)/3;
                       256*pt2^3 + 384*pt2^2*pt3 - 384*pt2^2 + 128*pt2*pt3^2 - 288*pt2*pt3 + 152*pt2 - 16*pt3^2 + 28*pt3 - 12,                                                                 28*pt2 - 32*pt2*pt3 + 128*pt2^2*pt3 - 144*pt2^2 + 128*pt2^3;
                       64*pt2*pt3 - (16*pt3)/3 - (224*pt2)/3 - 128*pt2^2*pt3 + 224*pt2^2 - (512*pt2^3)/3 + 16/3,                                                                             - (128*pt2^3)/3 + 32*pt2^2 - (16*pt2)/3;
                       384*pt2^2*pt3 + 512*pt2*pt3^2 - 448*pt2*pt3 + 128*pt3^3 - 224*pt3^2 + 96*pt3,                                                    128*pt2^3 + 512*pt2^2*pt3 - 224*pt2^2 + 384*pt2*pt3^2 - 448*pt2*pt3 + 96*pt2;
                       - 384*pt2^2*pt3 - 256*pt2*pt3^2 + 320*pt2*pt3 + 32*pt3^2 - 32*pt3,                                                                 64*pt2*pt3 - 32*pt2 - 256*pt2^2*pt3 + 160*pt2^2 - 128*pt2^3;
                       64*pt2*pt3 - 32*pt3 - 256*pt2*pt3^2 + 160*pt3^2 - 128*pt3^3,                                                             - 256*pt2^2*pt3 + 32*pt2^2 - 384*pt2*pt3^2 + 320*pt2*pt3 - 32*pt2];
        
        
        hess_matrix = [128*pt2^2 + 256*pt2*pt3 - 160*pt2 + 128*pt3^2 - 160*pt3 + 140/3,   128*pt2^2 + 256*pt2*pt3 - 160*pt2 + 128*pt3^2 - 160*pt3 + 140/3,   128*pt2^2 + 256*pt2*pt3 - 160*pt2 + 128*pt3^2 - 160*pt3 + 140/3,   128*pt2^2 + 256*pt2*pt3 - 160*pt2 + 128*pt3^2 - 160*pt3 + 140/3;
                       128*pt2^2 - 96*pt2 + 44/3,                                                     0,                                                     0,                                                     0;
                       0,                                                     0,                                                     0,                                 128*pt3^2 - 96*pt3 + 44/3;
                       256*pt2*pt3 - 64*pt3,                                 128*pt2^2 - 64*pt2 + 16/3,                                 128*pt2^2 - 64*pt2 + 16/3,                                                     0;
                       128*pt3^2 - 32*pt3,                             256*pt2*pt3 - 32*pt3 - 32*pt2 + 4,                             256*pt2*pt3 - 32*pt3 - 32*pt2 + 4,                                        128*pt2^2 - 32*pt2;
                       0,                                 128*pt3^2 - 64*pt3 + 16/3,                                 128*pt3^2 - 64*pt3 + 16/3,                                        256*pt2*pt3 - 64*pt2;
                       192*pt3 - 256*pt2*pt3 - 256*pt3^2, - 128*pt2^2 - 512*pt2*pt3 + 192*pt2 - 384*pt3^2 + 384*pt3 - 208/3, - 128*pt2^2 - 512*pt2*pt3 + 192*pt2 - 384*pt3^2 + 384*pt3 - 208/3, - 256*pt2^2 - 768*pt2*pt3 + 384*pt2 - 512*pt3^2 + 576*pt3 - 416/3;
                       128*pt3^2 - 32*pt3,                 256*pt2*pt3 - 288*pt3 - 32*pt2 + 384*pt3^2 + 28,                 256*pt2*pt3 - 288*pt3 - 32*pt2 + 384*pt3^2 + 28,     128*pt2^2 + 768*pt2*pt3 - 288*pt2 + 768*pt3^2 - 768*pt3 + 152;
                       0,                               - 128*pt3^2 + 64*pt3 - 16/3,                               - 128*pt3^2 + 64*pt3 - 16/3,              64*pt2 + 448*pt3 - 256*pt2*pt3 - 512*pt3^2 - 224/3;
                       - 512*pt2^2 - 768*pt2*pt3 + 576*pt2 - 256*pt3^2 + 384*pt3 - 416/3, - 384*pt2^2 - 512*pt2*pt3 + 384*pt2 - 128*pt3^2 + 192*pt3 - 208/3, - 384*pt2^2 - 512*pt2*pt3 + 384*pt2 - 128*pt3^2 + 192*pt3 - 208/3,                             192*pt2 - 256*pt2*pt3 - 256*pt2^2;
                       768*pt2^2 + 768*pt2*pt3 - 768*pt2 + 128*pt3^2 - 288*pt3 + 152,                 256*pt2*pt3 - 32*pt3 - 288*pt2 + 384*pt2^2 + 28,                 256*pt2*pt3 - 32*pt3 - 288*pt2 + 384*pt2^2 + 28,                                        128*pt2^2 - 32*pt2;
                       448*pt2 + 64*pt3 - 256*pt2*pt3 - 512*pt2^2 - 224/3,                               - 128*pt2^2 + 64*pt2 - 16/3,                               - 128*pt2^2 + 64*pt2 - 16/3,                                                     0;
                       768*pt2*pt3 - 448*pt3 + 512*pt3^2,     384*pt2^2 + 1024*pt2*pt3 - 448*pt2 + 384*pt3^2 - 448*pt3 + 96,     384*pt2^2 + 1024*pt2*pt3 - 448*pt2 + 384*pt3^2 - 448*pt3 + 96,                             768*pt2*pt3 - 448*pt2 + 512*pt2^2;
                       320*pt3 - 768*pt2*pt3 - 256*pt3^2,                 320*pt2 + 64*pt3 - 512*pt2*pt3 - 384*pt2^2 - 32,                 320*pt2 + 64*pt3 - 512*pt2*pt3 - 384*pt2^2 - 32,                                      - 256*pt2^2 + 64*pt2;
                       - 256*pt3^2 + 64*pt3,                 64*pt2 + 320*pt3 - 512*pt2*pt3 - 384*pt3^2 - 32,                 64*pt2 + 320*pt3 - 512*pt2*pt3 - 384*pt3^2 - 32,                             320*pt2 - 768*pt2*pt3 - 256*pt2^2];

         node_temp = [node(elem(i, 1), :); node(elem(i, 2), :); node(elem(i, 3), :)];
         node1 = 3/4*node(elem(i, 2), :) + 1/4*node(elem(i, 3), :);
         node1 = surfacedata.project(node1);
         node2 = 1/2*node(elem(i, 2), :) + 1/2*node(elem(i, 3), :);
         node2 = surfacedata.project(node2);
         node3 = 1/4*node(elem(i, 2), :) + 3/4*node(elem(i, 3), :);
         node3 = surfacedata.project(node3);
         node4 = 1/4*node(elem(i, 3), :) + 3/4*node(elem(i, 1), :);
         node4 = surfacedata.project(node4);
         node5 = 1/2*node(elem(i, 3), :) + 1/2*node(elem(i, 1), :);
         node5 = surfacedata.project(node5);
         node6 = 1/4*node(elem(i, 1), :) + 3/4*node(elem(i, 3), :);
         node6 = surfacedata.project(node6);
         node7 = 1/4*node(elem(i, 2), :) + 3/4*node(elem(i, 1), :);
         node7 = surfacedata.project(node7);
         node8 = 1/2*node(elem(i, 2), :) + 1/2*node(elem(i, 1), :);
         node8 = surfacedata.project(node8);
         node9 = 3/4*node(elem(i, 2), :) + 1/4*node(elem(i, 1), :);
         node9 = surfacedata.project(node9);
         node10 = 1/4*node(elem(i, 2), :) + 1/2*node(elem(i, 1), :) + 1/4*node(elem(i, 3), :);
         node10 = surfacedata.project(node10);
         node11 = 1/2*node(elem(i, 2), :) + 1/4*node(elem(i, 1), :) + 1/4*node(elem(i, 3), :);
         node11 = surfacedata.project(node11);
         node12 = 1/4*node(elem(i, 2), :) + 1/4*node(elem(i, 1), :) + 1/2*node(elem(i, 3), :);
         node12 = surfacedata.project(node12);
         node_temp(4:15, :) = [node1;
                             node2;
                             node3;
                             node4;
                             node5;
                             node6;
                             node7;
                             node8;
                             node9;
                             node10;
                             node11;
                             node12];
         for m=1:15
            %Derivatives of interpolation of geometry for building metric
            compx1 = compx1 + matrix(m, 1) * node_temp(m, 1);
            compy1 = compy1 + matrix(m, 1) * node_temp(m, 2);
            compz1 = compz1 + matrix(m, 1) * node_temp(m, 3);
            
            compx2 = compx2 + matrix(m, 2) * node_temp(m, 1);
            compy2 = compy2 + matrix(m, 2) * node_temp(m, 2);            
            compz2 = compz2 + matrix(m, 2) * node_temp(m, 3);

            hess_compx1 = hess_compx1 + hess_matrix(m, 1) * node_temp(m, 1);
            hess_compy1 = hess_compy1 + hess_matrix(m, 1) * node_temp(m, 2);
            hess_compz1 = hess_compz1 + hess_matrix(m, 1) * node_temp(m, 3);

            hess_compx2 = hess_compx2 + hess_matrix(m, 2) * node_temp(m, 1);
            hess_compy2 = hess_compy2 + hess_matrix(m, 2) * node_temp(m, 2);
            hess_compz2 = hess_compz2 + hess_matrix(m, 2) * node_temp(m, 3);

            hess_compx3 = hess_compx3 + hess_matrix(m, 3) * node_temp(m, 1);
            hess_compy3 = hess_compy3 + hess_matrix(m, 3) * node_temp(m, 2);
            hess_compz3 = hess_compz3 + hess_matrix(m, 3) * node_temp(m, 3);

            hess_compx4 = hess_compx4 + hess_matrix(m, 4) * node_temp(m, 1);
            hess_compy4 = hess_compy4 + hess_matrix(m, 4) * node_temp(m, 2);
            hess_compz4 = hess_compz4 + hess_matrix(m, 4) * node_temp(m, 3);

            point_comp1 = point_comp1 + point_matrix(m, 1) * node_temp(m, 1);
            point_comp2 = point_comp2 + point_matrix(m, 1) * node_temp(m, 2);
            point_comp3 = point_comp3 + point_matrix(m, 1) * node_temp(m, 3);
        end

        grad_ak_elem = [compx1, compx2; compy1 compy2; compz1 compz2]; 
        hess_ak_elem1 = [hess_compx1, hess_compx2; hess_compx3,hess_compx4];
        hess_ak_elem2 = [hess_compy1, hess_compy2; hess_compy3,hess_compy4];
        hess_ak_elem3 = [hess_compz1, hess_compz2; hess_compz3,hess_compz4];
        point_val = [point_comp1, point_comp2, point_comp3];
    case 5

        point_matrix = [-625*pt2^5/24 - 3125*pt2^4*pt3/24 + 625*pt2^4/8 - 3125*pt2^3*pt3^2/12 + 625*pt2^3*pt3/2 - 2125*pt2^3/24 - 3125*pt2^2*pt3^3/12 + 1875*pt2^2*pt3^2/4 - 2125*pt2^2*pt3/8 + 375*pt2^2/8 - 3125*pt2*pt3^4/24 + 625*pt2*pt3^3/2 - 2125*pt2*pt3^2/8 + 375*pt2*pt3/4 - 137*pt2/12 - 625*pt3^5/24 + 625*pt3^4/8 - 2125*pt3^3/24 + 375*pt3^2/8 - 137*pt3/12 + 1;
                        625*pt2^5/24 - 625*pt2^4/12 + 875*pt2^3/24 - 125*pt2^2/12 + pt2;
                        625*pt3^5/24 - 625*pt3^4/12 + 875*pt3^3/24 - 125*pt3^2/12 + pt3;
                        3125*pt2^4*pt3/24 - 625*pt2^3*pt3/4 + 1375*pt2^2*pt3/24 - 25*pt2*pt3/4;
                        3125*pt2^3*pt3^2/12 - 625*pt2^3*pt3/12 - 625*pt2^2*pt3^2/4 + 125*pt2^2*pt3/4 + 125*pt2*pt3^2/6 - 25*pt2*pt3/6;
                        3125*pt2^2*pt3^3/12 - 625*pt2^2*pt3^2/4 + 125*pt2^2*pt3/6 - 625*pt2*pt3^3/12 + 125*pt2*pt3^2/4 - 25*pt2*pt3/6;
                        3125*pt2*pt3^4/24 - 625*pt2*pt3^3/4 + 1375*pt2*pt3^2/24 - 25*pt2*pt3/4;
                        3125*pt2^4*pt3/24 + 3125*pt2^3*pt3^2/6 - 4375*pt2^3*pt3/12 + 3125*pt2^2*pt3^3/4 - 4375*pt2^2*pt3^2/4 + 8875*pt2^2*pt3/24 + 3125*pt2*pt3^4/6 - 4375*pt2*pt3^3/4 + 8875*pt2*pt3^2/12 - 1925*pt2*pt3/12 + 3125*pt3^5/24 - 4375*pt3^4/12 + 8875*pt3^3/24 - 1925*pt3^2/12 + 25*pt3;
                        -3125*pt2^3*pt3^2/12 + 625*pt2^3*pt3/12 - 3125*pt2^2*pt3^3/4 + 3125*pt2^2*pt3^2/4 - 125*pt2^2*pt3 - 3125*pt2*pt3^4/4 + 5625*pt2*pt3^3/4 - 8875*pt2*pt3^2/12 + 1175*pt2*pt3/12 - 3125*pt3^5/12 + 8125*pt3^4/12 - 7375*pt3^3/12 + 2675*pt3^2/12 - 25*pt3;
                        3125*pt2^2*pt3^3/12 - 625*pt2^2*pt3^2/4 + 125*pt2^2*pt3/6 + 3125*pt2*pt3^4/6 - 3125*pt2*pt3^3/4 + 3875*pt2*pt3^2/12 - 75*pt2*pt3/2 + 3125*pt3^5/12 - 625*pt3^4 + 6125*pt3^3/12 - 325*pt3^2/2 + 50*pt3/3;
                        -3125*pt2*pt3^4/24 + 625*pt2*pt3^3/4 - 1375*pt2*pt3^2/24 + 25*pt2*pt3/4 - 3125*pt3^5/24 + 6875*pt3^4/24 - 5125*pt3^3/24 + 1525*pt3^2/24 - 25*pt3/4;
                        3125*pt2^5/24 + 3125*pt2^4*pt3/6 - 4375*pt2^4/12 + 3125*pt2^3*pt3^2/4 - 4375*pt2^3*pt3/4 + 8875*pt2^3/24 + 3125*pt2^2*pt3^3/6 - 4375*pt2^2*pt3^2/4 + 8875*pt2^2*pt3/12 - 1925*pt2^2/12 + 3125*pt2*pt3^4/24 - 4375*pt2*pt3^3/12 + 8875*pt2*pt3^2/24 - 1925*pt2*pt3/12 + 25*pt2;
                        -3125*pt2^5/12 - 3125*pt2^4*pt3/4 + 8125*pt2^4/12 - 3125*pt2^3*pt3^2/4 + 5625*pt2^3*pt3/4 - 7375*pt2^3/12 - 3125*pt2^2*pt3^3/12 + 3125*pt2^2*pt3^2/4 - 8875*pt2^2*pt3/12 + 2675*pt2^2/12 + 625*pt2*pt3^3/12 - 125*pt2*pt3^2 + 1175*pt2*pt3/12 - 25*pt2;
                        3125*pt2^5/12 + 3125*pt2^4*pt3/6 - 625*pt2^4 + 3125*pt2^3*pt3^2/12 - 3125*pt2^3*pt3/4 + 6125*pt2^3/12 - 625*pt2^2*pt3^2/4 + 3875*pt2^2*pt3/12 - 325*pt2^2/2 + 125*pt2*pt3^2/6 - 75*pt2*pt3/2 + 50*pt2/3;
                        -3125*pt2^5/24 - 3125*pt2^4*pt3/24 + 6875*pt2^4/24 + 625*pt2^3*pt3/4 - 5125*pt2^3/24 - 1375*pt2^2*pt3/24 + 1525*pt2^2/24 + 25*pt2*pt3/4 - 25*pt2/4;
                        -3125*pt2^4*pt3/6 - 3125*pt2^3*pt3^2/2 + 1250*pt2^3*pt3 - 3125*pt2^2*pt3^3/2 + 2500*pt2^2*pt3^2 - 5875*pt2^2*pt3/6 - 3125*pt2*pt3^4/6 + 1250*pt2*pt3^3 - 5875*pt2*pt3^2/6 + 250*pt2*pt3;
                        3125*pt2^4*pt3/4 + 3125*pt2^3*pt3^2/2 - 3125*pt2^3*pt3/2 + 3125*pt2^2*pt3^3/4 - 6875*pt2^2*pt3^2/4 + 3625*pt2^2*pt3/4 - 625*pt2*pt3^3/4 + 1125*pt2*pt3^2/4 - 125*pt2*pt3;
                        -3125*pt2^4*pt3/6 - 3125*pt2^3*pt3^2/6 + 2500*pt2^3*pt3/3 + 625*pt2^2*pt3^2/2 - 2125*pt2^2*pt3/6 - 125*pt2*pt3^2/3 + 125*pt2*pt3/3;
                        3125*pt2^3*pt3^2/4 - 625*pt2^3*pt3/4 + 3125*pt2^2*pt3^3/2 - 6875*pt2^2*pt3^2/4 + 1125*pt2^2*pt3/4 + 3125*pt2*pt3^4/4 - 3125*pt2*pt3^3/2 + 3625*pt2*pt3^2/4 - 125*pt2*pt3;
                        -3125*pt2^3*pt3^2/4 + 625*pt2^3*pt3/4 - 3125*pt2^2*pt3^3/4 + 4375*pt2^2*pt3^2/4 - 375*pt2^2*pt3/2 + 625*pt2*pt3^3/4 - 375*pt2*pt3^2/2 + 125*pt2*pt3/4;
                        -3125*pt2^2*pt3^3/6 + 625*pt2^2*pt3^2/2 - 125*pt2^2*pt3/3 - 3125*pt2*pt3^4/6 + 2500*pt2*pt3^3/3 - 2125*pt2*pt3^2/6 + 125*pt2*pt3/3]; 

            matrix = [   - (3125*pt2^4)/24 - (3125*pt2^3*pt3)/6 + (625*pt2^3)/2 - (3125*pt2^2*pt3^2)/4 + (1875*pt2^2*pt3)/2 - (2125*pt2^2)/8 - (3125*pt2*pt3^3)/6 + (1875*pt2*pt3^2)/2 - (2125*pt2*pt3)/4 + (375*pt2)/4 - (3125*pt3^4)/24 + (625*pt3^3)/2 - (2125*pt3^2)/8 + (375*pt3)/4 - 137/12,    - (3125*pt2^4)/24 - (3125*pt2^3*pt3)/6 + (625*pt2^3)/2 - (3125*pt2^2*pt3^2)/4 + (1875*pt2^2*pt3)/2 - (2125*pt2^2)/8 - (3125*pt2*pt3^3)/6 + (1875*pt2*pt3^2)/2 - (2125*pt2*pt3)/4 + (375*pt2)/4 - (3125*pt3^4)/24 + (625*pt3^3)/2 - (2125*pt3^2)/8 + (375*pt3)/4 - 137/12;
                                                                                                                                                                          (3125*pt2^4)/24 - (625*pt2^3)/3 + (875*pt2^2)/8 - (125*pt2)/6 + 1,                                                                                                                                                                                                                                   0;
                                                                                                                                                                                                                                  0,                                                                                                                                                                           (3125*pt3^4)/24 - (625*pt3^3)/3 + (875*pt3^2)/8 - (125*pt3)/6 + 1;
                                                                                                                                                                        (3125*pt3*pt2^3)/6 - (1875*pt3*pt2^2)/4 + (1375*pt3*pt2)/12 - (25*pt3)/4,                                                                                                                                                                              (3125*pt2^4)/24 - (625*pt2^3)/4 + (1375*pt2^2)/24 - (25*pt2)/4;
                                                                                                                                            (3125*pt2^2*pt3^2)/4 - (625*pt2^2*pt3)/4 - (625*pt2*pt3^2)/2 + (125*pt2*pt3)/2 + (125*pt3^2)/6 - (25*pt3)/6,                                                                                                                                                (125*pt2*pt3)/3 - (25*pt2)/6 - (625*pt2^2*pt3)/2 + (3125*pt2^3*pt3)/6 + (125*pt2^2)/4 - (625*pt2^3)/12;
                                                                                                                                               (125*pt2*pt3)/3 - (25*pt3)/6 - (625*pt2*pt3^2)/2 + (3125*pt2*pt3^3)/6 + (125*pt3^2)/4 - (625*pt3^3)/12,                                                                                                                                             (3125*pt2^2*pt3^2)/4 - (625*pt2^2*pt3)/2 + (125*pt2^2)/6 - (625*pt2*pt3^2)/4 + (125*pt2*pt3)/2 - (25*pt2)/6;
                                                                                                                                                                            (3125*pt3^4)/24 - (625*pt3^3)/4 + (1375*pt3^2)/24 - (25*pt3)/4,                                                                                                                                                                          (3125*pt2*pt3^3)/6 - (1875*pt2*pt3^2)/4 + (1375*pt2*pt3)/12 - (25*pt2)/4;
                                                                  (3125*pt2^3*pt3)/6 + (3125*pt2^2*pt3^2)/2 - (4375*pt2^2*pt3)/4 + (3125*pt2*pt3^3)/2 - (4375*pt2*pt3^2)/2 + (8875*pt2*pt3)/12 + (3125*pt3^4)/6 - (4375*pt3^3)/4 + (8875*pt3^2)/12 - (1925*pt3)/12, (3125*pt2^4)/24 + (3125*pt2^3*pt3)/3 - (4375*pt2^3)/12 + (9375*pt2^2*pt3^2)/4 - (4375*pt2^2*pt3)/2 + (8875*pt2^2)/24 + (6250*pt2*pt3^3)/3 - (13125*pt2*pt3^2)/4 + (8875*pt2*pt3)/6 - (1925*pt2)/12 + (15625*pt3^4)/24 - (4375*pt3^3)/3 + (8875*pt3^2)/8 - (1925*pt3)/6 + 25;
                                                                                         - (3125*pt2^2*pt3^2)/4 + (625*pt2^2*pt3)/4 - (3125*pt2*pt3^3)/2 + (3125*pt2*pt3^2)/2 - 250*pt2*pt3 - (3125*pt3^4)/4 + (5625*pt3^3)/4 - (8875*pt3^2)/12 + (1175*pt3)/12,                          - (3125*pt2^3*pt3)/6 + (625*pt2^3)/12 - (9375*pt2^2*pt3^2)/4 + (3125*pt2^2*pt3)/2 - 125*pt2^2 - 3125*pt2*pt3^3 + (16875*pt2*pt3^2)/4 - (8875*pt2*pt3)/6 + (1175*pt2)/12 - (15625*pt3^4)/12 + (8125*pt3^3)/3 - (7375*pt3^2)/4 + (2675*pt3)/6 - 25;
                                                                                                                              (125*pt2*pt3)/3 - (75*pt3)/2 - (625*pt2*pt3^2)/2 + (3125*pt2*pt3^3)/6 + (3875*pt3^2)/12 - (3125*pt3^3)/4 + (3125*pt3^4)/6,                                                                (3125*pt2^2*pt3^2)/4 - (625*pt2^2*pt3)/2 + (125*pt2^2)/6 + (6250*pt2*pt3^3)/3 - (9375*pt2*pt3^2)/4 + (3875*pt2*pt3)/6 - (75*pt2)/2 + (15625*pt3^4)/12 - 2500*pt3^3 + (6125*pt3^2)/4 - 325*pt3 + 50/3;
                                                                                                                                                                           - (3125*pt3^4)/24 + (625*pt3^3)/4 - (1375*pt3^2)/24 + (25*pt3)/4,                                                                                                      (25*pt2)/4 + (1525*pt3)/12 - (1375*pt2*pt3)/12 + (1875*pt2*pt3^2)/4 - (3125*pt2*pt3^3)/6 - (5125*pt3^2)/8 + (6875*pt3^3)/6 - (15625*pt3^4)/24 - 25/4;
                                    (15625*pt2^4)/24 + (6250*pt2^3*pt3)/3 - (4375*pt2^3)/3 + (9375*pt2^2*pt3^2)/4 - (13125*pt2^2*pt3)/4 + (8875*pt2^2)/8 + (3125*pt2*pt3^3)/3 - (4375*pt2*pt3^2)/2 + (8875*pt2*pt3)/6 - (1925*pt2)/6 + (3125*pt3^4)/24 - (4375*pt3^3)/12 + (8875*pt3^2)/24 - (1925*pt3)/12 + 25,                                                                    (3125*pt2^4)/6 + (3125*pt2^3*pt3)/2 - (4375*pt2^3)/4 + (3125*pt2^2*pt3^2)/2 - (4375*pt2^2*pt3)/2 + (8875*pt2^2)/12 + (3125*pt2*pt3^3)/6 - (4375*pt2*pt3^2)/4 + (8875*pt2*pt3)/12 - (1925*pt2)/12;
                                     - (15625*pt2^4)/12 - 3125*pt2^3*pt3 + (8125*pt2^3)/3 - (9375*pt2^2*pt3^2)/4 + (16875*pt2^2*pt3)/4 - (7375*pt2^2)/4 - (3125*pt2*pt3^3)/6 + (3125*pt2*pt3^2)/2 - (8875*pt2*pt3)/6 + (2675*pt2)/6 + (625*pt3^3)/12 - 125*pt3^2 + (1175*pt3)/12 - 25,                                                                                          - (3125*pt2^4)/4 - (3125*pt2^3*pt3)/2 + (5625*pt2^3)/4 - (3125*pt2^2*pt3^2)/4 + (3125*pt2^2*pt3)/2 - (8875*pt2^2)/12 + (625*pt2*pt3^2)/4 - 250*pt2*pt3 + (1175*pt2)/12;
                                                               (15625*pt2^4)/12 + (6250*pt2^3*pt3)/3 - 2500*pt2^3 + (3125*pt2^2*pt3^2)/4 - (9375*pt2^2*pt3)/4 + (6125*pt2^2)/4 - (625*pt2*pt3^2)/2 + (3875*pt2*pt3)/6 - 325*pt2 + (125*pt3^2)/6 - (75*pt3)/2 + 50/3,                                                                                                                               (125*pt2*pt3)/3 - (75*pt2)/2 - (625*pt2^2*pt3)/2 + (3125*pt2^3*pt3)/6 + (3875*pt2^2)/12 - (3125*pt2^3)/4 + (3125*pt2^4)/6;
                                                                                                     (1525*pt2)/12 + (25*pt3)/4 - (1375*pt2*pt3)/12 + (1875*pt2^2*pt3)/4 - (3125*pt2^3*pt3)/6 - (5125*pt2^2)/8 + (6875*pt2^3)/6 - (15625*pt2^4)/24 - 25/4,                                                                                                                                                                            - (3125*pt2^4)/24 + (625*pt2^3)/4 - (1375*pt2^2)/24 + (25*pt2)/4;
                                                                                         - (6250*pt2^3*pt3)/3 - (9375*pt2^2*pt3^2)/2 + 3750*pt2^2*pt3 - 3125*pt2*pt3^3 + 5000*pt2*pt3^2 - (5875*pt2*pt3)/3 - (3125*pt3^4)/6 + 1250*pt3^3 - (5875*pt3^2)/6 + 250*pt3,                                                                                          - (3125*pt2^4)/6 - 3125*pt2^3*pt3 + 1250*pt2^3 - (9375*pt2^2*pt3^2)/2 + 5000*pt2^2*pt3 - (5875*pt2^2)/6 - (6250*pt2*pt3^3)/3 + 3750*pt2*pt3^2 - (5875*pt2*pt3)/3 + 250*pt2;
                                                                                               3125*pt2^3*pt3 + (9375*pt2^2*pt3^2)/2 - (9375*pt2^2*pt3)/2 + (3125*pt2*pt3^3)/2 - (6875*pt2*pt3^2)/2 + (3625*pt2*pt3)/2 - (625*pt3^3)/4 + (1125*pt3^2)/4 - 125*pt3,                                                                                                 (3125*pt2^4)/4 + 3125*pt2^3*pt3 - (3125*pt2^3)/2 + (9375*pt2^2*pt3^2)/4 - (6875*pt2^2*pt3)/2 + (3625*pt2^2)/4 - (1875*pt2*pt3^2)/4 + (1125*pt2*pt3)/2 - 125*pt2;
                                                                                                                              - (6250*pt2^3*pt3)/3 - (3125*pt2^2*pt3^2)/2 + 2500*pt2^2*pt3 + 625*pt2*pt3^2 - (2125*pt2*pt3)/3 - (125*pt3^2)/3 + (125*pt3)/3,                                                                                                                                   (125*pt2)/3 - (250*pt2*pt3)/3 + 625*pt2^2*pt3 - (3125*pt2^3*pt3)/3 - (2125*pt2^2)/6 + (2500*pt2^3)/3 - (3125*pt2^4)/6;
                                                                                                (9375*pt2^2*pt3^2)/4 - (1875*pt2^2*pt3)/4 + 3125*pt2*pt3^3 - (6875*pt2*pt3^2)/2 + (1125*pt2*pt3)/2 + (3125*pt3^4)/4 - (3125*pt3^3)/2 + (3625*pt3^2)/4 - 125*pt3,                                                                                                (3125*pt2^3*pt3)/2 - (625*pt2^3)/4 + (9375*pt2^2*pt3^2)/2 - (6875*pt2^2*pt3)/2 + (1125*pt2^2)/4 + 3125*pt2*pt3^3 - (9375*pt2*pt3^2)/2 + (3625*pt2*pt3)/2 - 125*pt2;
                                                                                                            - (9375*pt2^2*pt3^2)/4 + (1875*pt2^2*pt3)/4 - (3125*pt2*pt3^3)/2 + (4375*pt2*pt3^2)/2 - 375*pt2*pt3 + (625*pt3^3)/4 - (375*pt3^2)/2 + (125*pt3)/4,                                                                                                             - (3125*pt2^3*pt3)/2 + (625*pt2^3)/4 - (9375*pt2^2*pt3^2)/4 + (4375*pt2^2*pt3)/2 - (375*pt2^2)/2 + (1875*pt2*pt3^2)/4 - 375*pt2*pt3 + (125*pt2)/4;
                                      (125*pt3)/3 - (250*pt2*pt3)/3 + 625*pt2*pt3^2 - (3125*pt2*pt3^3)/3 - (2125*pt3^2)/6 + (2500*pt3^3)/3 - (3125*pt3^4)/6,                                                                                                                               - (3125*pt2^2*pt3^2)/2 + 625*pt2^2*pt3 - (125*pt2^2)/3 - (6250*pt2*pt3^3)/3 + 2500*pt2*pt3^2 - (2125*pt2*pt3)/3 + (125*pt2)/3];

                            
                 hess_matrix = [- (3125*pt2^3)/6 - (3125*pt2^2*pt3)/2 + (1875*pt2^2)/2 - (3125*pt2*pt3^2)/2 + 1875*pt2*pt3 - (2125*pt2)/4 - (3125*pt3^3)/6 + (1875*pt3^2)/2 - (2125*pt3)/4 + 375/4, - (3125*pt2^3)/6 - (3125*pt2^2*pt3)/2 + (1875*pt2^2)/2 - (3125*pt2*pt3^2)/2 + 1875*pt2*pt3 - (2125*pt2)/4 - (3125*pt3^3)/6 + (1875*pt3^2)/2 - (2125*pt3)/4 + 375/4, - (3125*pt2^3)/6 - (3125*pt2^2*pt3)/2 + (1875*pt2^2)/2 - (3125*pt2*pt3^2)/2 + 1875*pt2*pt3 - (2125*pt2)/4 - (3125*pt3^3)/6 + (1875*pt3^2)/2 - (2125*pt3)/4 + 375/4, - (3125*pt2^3)/6 - (3125*pt2^2*pt3)/2 + (1875*pt2^2)/2 - (3125*pt2*pt3^2)/2 + 1875*pt2*pt3 - (2125*pt2)/4 - (3125*pt3^3)/6 + (1875*pt3^2)/2 - (2125*pt3)/4 + 375/4;
                                                                                                (3125*pt2^3)/6 - 625*pt2^2 + (875*pt2)/4 - 125/6,                                                                                                                                          0,                                                                                                                                          0,                                                                                                                                          0;
                                                                                                                                         0,                                                                                                                                          0,                                                                                                                                          0,                                                                                                 (3125*pt3^3)/6 - 625*pt3^2 + (875*pt3)/4 - 125/6;
                                                                                               (3125*pt3*pt2^2)/2 - (1875*pt3*pt2)/2 + (1375*pt3)/12,                                                                                           (3125*pt2^3)/6 - (1875*pt2^2)/4 + (1375*pt2)/12 - 25/4,                                                                                           (3125*pt2^3)/6 - (1875*pt2^2)/4 + (1375*pt2)/12 - 25/4,                                                                                                                                          0;
                                                                                    (125*pt3)/2 - (625*pt2*pt3)/2 + (3125*pt2*pt3^2)/2 - (625*pt3^2)/2,                                                                      (125*pt2)/2 + (125*pt3)/3 - 625*pt2*pt3 + (3125*pt2^2*pt3)/2 - (625*pt2^2)/4 - 25/6,                                                                      (125*pt2)/2 + (125*pt3)/3 - 625*pt2*pt3 + (3125*pt2^2*pt3)/2 - (625*pt2^2)/4 - 25/6,                                                                                                     (3125*pt2^3)/6 - (625*pt2^2)/2 + (125*pt2)/3;
                                                                                                    (3125*pt3^3)/6 - (625*pt3^2)/2 + (125*pt3)/3,                                                                      (125*pt2)/3 + (125*pt3)/2 - 625*pt2*pt3 + (3125*pt2*pt3^2)/2 - (625*pt3^2)/4 - 25/6,                                                                      (125*pt2)/3 + (125*pt3)/2 - 625*pt2*pt3 + (3125*pt2*pt3^2)/2 - (625*pt3^2)/4 - 25/6,                                                                                     (125*pt2)/2 - (625*pt2*pt3)/2 + (3125*pt2^2*pt3)/2 - (625*pt2^2)/2;
                                                                                                                                         0,                                                                                           (3125*pt3^3)/6 - (1875*pt3^2)/4 + (1375*pt3)/12 - 25/4,                                                                                           (3125*pt3^3)/6 - (1875*pt3^2)/4 + (1375*pt3)/12 - 25/4,                                                                                                (3125*pt2*pt3^2)/2 - (1875*pt2*pt3)/2 + (1375*pt2)/12;
                                                    (3125*pt2^2*pt3)/2 + 3125*pt2*pt3^2 - (4375*pt2*pt3)/2 + (3125*pt3^3)/2 - (4375*pt3^2)/2 + (8875*pt3)/12,   (3125*pt2^3)/6 + 3125*pt2^2*pt3 - (4375*pt2^2)/4 + (9375*pt2*pt3^2)/2 - 4375*pt2*pt3 + (8875*pt2)/12 + (6250*pt3^3)/3 - (13125*pt3^2)/4 + (8875*pt3)/6 - 1925/12,   (3125*pt2^3)/6 + 3125*pt2^2*pt3 - (4375*pt2^2)/4 + (9375*pt2*pt3^2)/2 - 4375*pt2*pt3 + (8875*pt2)/12 + (6250*pt3^3)/3 - (13125*pt3^2)/4 + (8875*pt3)/6 - 1925/12,    (3125*pt2^3)/3 + (9375*pt2^2*pt3)/2 - (4375*pt2^2)/2 + 6250*pt2*pt3^2 - (13125*pt2*pt3)/2 + (8875*pt2)/6 + (15625*pt3^3)/6 - 4375*pt3^2 + (8875*pt3)/4 - 1925/6;
                                                                        (625*pt2*pt3)/2 - 250*pt3 - (3125*pt2*pt3^2)/2 + (3125*pt3^2)/2 - (3125*pt3^3)/2,                       - (3125*pt2^2*pt3)/2 + (625*pt2^2)/4 - (9375*pt2*pt3^2)/2 + 3125*pt2*pt3 - 250*pt2 - 3125*pt3^3 + (16875*pt3^2)/4 - (8875*pt3)/6 + 1175/12,                       - (3125*pt2^2*pt3)/2 + (625*pt2^2)/4 - (9375*pt2*pt3^2)/2 + 3125*pt2*pt3 - 250*pt2 - 3125*pt3^3 + (16875*pt3^2)/4 - (8875*pt3)/6 + 1175/12,  - (3125*pt2^3)/6 - (9375*pt2^2*pt3)/2 + (3125*pt2^2)/2 - 9375*pt2*pt3^2 + (16875*pt2*pt3)/2 - (8875*pt2)/6 - (15625*pt3^3)/3 + 8125*pt3^2 - (7375*pt3)/2 + 2675/6;
                                                                                                    (3125*pt3^3)/6 - (625*pt3^2)/2 + (125*pt3)/3,                                                     (125*pt2)/3 + (3875*pt3)/6 - 625*pt2*pt3 + (3125*pt2*pt3^2)/2 - (9375*pt3^2)/4 + (6250*pt3^3)/3 - 75/2,                                                     (125*pt2)/3 + (3875*pt3)/6 - 625*pt2*pt3 + (3125*pt2*pt3^2)/2 - (9375*pt3^2)/4 + (6250*pt3^3)/3 - 75/2,                        (3125*pt2^2*pt3)/2 - (625*pt2^2)/2 + 6250*pt2*pt3^2 - (9375*pt2*pt3)/2 + (3875*pt2)/6 + (15625*pt3^3)/3 - 7500*pt3^2 + (6125*pt3)/2 - 325;
                                                                                                                                         0,                                                                                         - (3125*pt3^3)/6 + (1875*pt3^2)/4 - (1375*pt3)/12 + 25/4,                                                                                         - (3125*pt3^3)/6 + (1875*pt3^2)/4 - (1375*pt3)/12 + 25/4,                                          (1875*pt2*pt3)/2 - (5125*pt3)/4 - (1375*pt2)/12 - (3125*pt2*pt3^2)/2 + (6875*pt3^2)/2 - (15625*pt3^3)/6 + 1525/12;
                             (15625*pt2^3)/6 + 6250*pt2^2*pt3 - 4375*pt2^2 + (9375*pt2*pt3^2)/2 - (13125*pt2*pt3)/2 + (8875*pt2)/4 + (3125*pt3^3)/3 - (4375*pt3^2)/2 + (8875*pt3)/6 - 1925/6,   (6250*pt2^3)/3 + (9375*pt2^2*pt3)/2 - (13125*pt2^2)/4 + 3125*pt2*pt3^2 - 4375*pt2*pt3 + (8875*pt2)/6 + (3125*pt3^3)/6 - (4375*pt3^2)/4 + (8875*pt3)/12 - 1925/12,   (6250*pt2^3)/3 + (9375*pt2^2*pt3)/2 - (13125*pt2^2)/4 + 3125*pt2*pt3^2 - 4375*pt2*pt3 + (8875*pt2)/6 + (3125*pt3^3)/6 - (4375*pt3^2)/4 + (8875*pt3)/12 - 1925/12,                                                     (3125*pt2^3)/2 + 3125*pt2^2*pt3 - (4375*pt2^2)/2 + (3125*pt2*pt3^2)/2 - (4375*pt2*pt3)/2 + (8875*pt2)/12;
                            - (15625*pt2^3)/3 - 9375*pt2^2*pt3 + 8125*pt2^2 - (9375*pt2*pt3^2)/2 + (16875*pt2*pt3)/2 - (7375*pt2)/2 - (3125*pt3^3)/6 + (3125*pt3^2)/2 - (8875*pt3)/6 + 2675/6,                       - 3125*pt2^3 - (9375*pt2^2*pt3)/2 + (16875*pt2^2)/4 - (3125*pt2*pt3^2)/2 + 3125*pt2*pt3 - (8875*pt2)/6 + (625*pt3^2)/4 - 250*pt3 + 1175/12,                       - 3125*pt2^3 - (9375*pt2^2*pt3)/2 + (16875*pt2^2)/4 - (3125*pt2*pt3^2)/2 + 3125*pt2*pt3 - (8875*pt2)/6 + (625*pt3^2)/4 - 250*pt3 + 1175/12,                                                                         (625*pt2*pt3)/2 - 250*pt2 - (3125*pt2^2*pt3)/2 + (3125*pt2^2)/2 - (3125*pt2^3)/2;
                             (15625*pt2^3)/3 + 6250*pt2^2*pt3 - 7500*pt2^2 + (3125*pt2*pt3^2)/2 - (9375*pt2*pt3)/2 + (6125*pt2)/2 - (625*pt3^2)/2 + (3875*pt3)/6 - 325,                                                     (3875*pt2)/6 + (125*pt3)/3 - 625*pt2*pt3 + (3125*pt2^2*pt3)/2 - (9375*pt2^2)/4 + (6250*pt2^3)/3 - 75/2,                                                     (3875*pt2)/6 + (125*pt3)/3 - 625*pt2*pt3 + (3125*pt2^2*pt3)/2 - (9375*pt2^2)/4 + (6250*pt2^3)/3 - 75/2,                                                                                                     (3125*pt2^3)/6 - (625*pt2^2)/2 + (125*pt2)/3;
                                         (1875*pt2*pt3)/2 - (1375*pt3)/12 - (5125*pt2)/4 - (3125*pt2^2*pt3)/2 + (6875*pt2^2)/2 - (15625*pt2^3)/6 + 1525/12,                                                                                         - (3125*pt2^3)/6 + (1875*pt2^2)/4 - (1375*pt2)/12 + 25/4,                                                                                         - (3125*pt2^3)/6 + (1875*pt2^2)/4 - (1375*pt2)/12 + 25/4,                                                                                                                                          0;
                                                                   - 6250*pt2^2*pt3 - 9375*pt2*pt3^2 + 7500*pt2*pt3 - 3125*pt3^3 + 5000*pt3^2 - (5875*pt3)/3,                  - (6250*pt2^3)/3 - 9375*pt2^2*pt3 + 3750*pt2^2 - 9375*pt2*pt3^2 + 10000*pt2*pt3 - (5875*pt2)/3 - (6250*pt3^3)/3 + 3750*pt3^2 - (5875*pt3)/3 + 250,                  - (6250*pt2^3)/3 - 9375*pt2^2*pt3 + 3750*pt2^2 - 9375*pt2*pt3^2 + 10000*pt2*pt3 - (5875*pt2)/3 - (6250*pt3^3)/3 + 3750*pt3^2 - (5875*pt3)/3 + 250,                                                                    - 3125*pt2^3 - 9375*pt2^2*pt3 + 5000*pt2^2 - 6250*pt2*pt3^2 + 7500*pt2*pt3 - (5875*pt2)/3;
                                                             9375*pt2^2*pt3 + 9375*pt2*pt3^2 - 9375*pt2*pt3 + (3125*pt3^3)/2 - (6875*pt3^2)/2 + (3625*pt3)/2,                            3125*pt2^3 + 9375*pt2^2*pt3 - (9375*pt2^2)/2 + (9375*pt2*pt3^2)/2 - 6875*pt2*pt3 + (3625*pt2)/2 - (1875*pt3^2)/4 + (1125*pt3)/2 - 125,                            3125*pt2^3 + 9375*pt2^2*pt3 - (9375*pt2^2)/2 + (9375*pt2*pt3^2)/2 - 6875*pt2*pt3 + (3625*pt2)/2 - (1875*pt3^2)/4 + (1125*pt3)/2 - 125,                                                                       (1125*pt2)/2 - (1875*pt2*pt3)/2 + (9375*pt2^2*pt3)/2 - (6875*pt2^2)/2 + 3125*pt2^3;
                                                                               - 6250*pt2^2*pt3 - 3125*pt2*pt3^2 + 5000*pt2*pt3 + 625*pt3^2 - (2125*pt3)/3,                                                           1250*pt2*pt3 - (250*pt3)/3 - (2125*pt2)/3 - 3125*pt2^2*pt3 + 2500*pt2^2 - (6250*pt2^3)/3 + 125/3,                                                           1250*pt2*pt3 - (250*pt3)/3 - (2125*pt2)/3 - 3125*pt2^2*pt3 + 2500*pt2^2 - (6250*pt2^3)/3 + 125/3,                                                                                                       - (3125*pt2^3)/3 + 625*pt2^2 - (250*pt2)/3;
                                                                      (1125*pt3)/2 - (1875*pt2*pt3)/2 + (9375*pt2*pt3^2)/2 - (6875*pt3^2)/2 + 3125*pt3^3,                            (9375*pt2^2*pt3)/2 - (1875*pt2^2)/4 + 9375*pt2*pt3^2 - 6875*pt2*pt3 + (1125*pt2)/2 + 3125*pt3^3 - (9375*pt3^2)/2 + (3625*pt3)/2 - 125,                            (9375*pt2^2*pt3)/2 - (1875*pt2^2)/4 + 9375*pt2*pt3^2 - 6875*pt2*pt3 + (1125*pt2)/2 + 3125*pt3^3 - (9375*pt3^2)/2 + (3625*pt3)/2 - 125,                                                              (3125*pt2^3)/2 + 9375*pt2^2*pt3 - (6875*pt2^2)/2 + 9375*pt2*pt3^2 - 9375*pt2*pt3 + (3625*pt2)/2;
                                                                       (1875*pt2*pt3)/2 - 375*pt3 - (9375*pt2*pt3^2)/2 + (4375*pt3^2)/2 - (3125*pt3^3)/2,                                         - (9375*pt2^2*pt3)/2 + (1875*pt2^2)/4 - (9375*pt2*pt3^2)/2 + 4375*pt2*pt3 - 375*pt2 + (1875*pt3^2)/4 - 375*pt3 + 125/4,                                         - (9375*pt2^2*pt3)/2 + (1875*pt2^2)/4 - (9375*pt2*pt3^2)/2 + 4375*pt2*pt3 - 375*pt2 + (1875*pt3^2)/4 - 375*pt3 + 125/4,                                                                        (1875*pt2*pt3)/2 - 375*pt2 - (9375*pt2^2*pt3)/2 + (4375*pt2^2)/2 - (3125*pt2^3)/2;
                                                                                                      - (3125*pt3^3)/3 + 625*pt3^2 - (250*pt3)/3,                                                           1250*pt2*pt3 - (2125*pt3)/3 - (250*pt2)/3 - 3125*pt2*pt3^2 + 2500*pt3^2 - (6250*pt3^3)/3 + 125/3,                                                           1250*pt2*pt3 - (2125*pt3)/3 - (250*pt2)/3 - 3125*pt2*pt3^2 + 2500*pt3^2 - (6250*pt3^3)/3 + 125/3,                                                                                - 3125*pt2^2*pt3 + 625*pt2^2 - 6250*pt2*pt3^2 + 5000*pt2*pt3 - (2125*pt2)/3];
                 node_temp = [node(elem(i, 1), :); node(elem(i, 2), :); node(elem(i, 3), :)];
                 node1 = 4/5*node(elem(i, 2), :) + 1/5*node(elem(i, 3), :);
                 node1 = node1 /norm(node1, 2);
                 node2 = 3/5*node(elem(i, 2), :) + 2/5*node(elem(i, 3), :);
                 node2 = node2 /norm(node2, 2);
                 node3 = 2/5*node(elem(i, 2), :) + 3/5*node(elem(i, 3), :);
                 node3 = node3 /norm(node3, 2);
                 node4 = 1/5*node(elem(i, 2), :) + 4/5*node(elem(i, 3), :);
                 node4 = node4 /norm(node4, 2); 
                 
                 
                 node5 = 1/5*node(elem(i, 1), :) + 4/5*node(elem(i, 3), :);
                 node5 = node5 /norm(node5, 2);
                 node6 = 2/5*node(elem(i, 1), :) + 3/5*node(elem(i, 3), :);
                 node6 = node6 /norm(node6, 2);
                 node7 = 3/5*node(elem(i, 1), :) + 2/5*node(elem(i, 3), :);
                 node7 = node7 /norm(node7, 2);
                 node8 = 4/5*node(elem(i, 1), :) + 1/5*node(elem(i, 3), :);
                 node8 = node8 /norm(node8, 2);


                 node9 = 1/5*node(elem(i, 2), :) + 4/5*node(elem(i, 1), :);
                 node9 = node9 /norm(node9, 2);
                 node10 = 2/5*node(elem(i, 2), :) + 3/5*node(elem(i, 1), :);
                 node10 = node10 /norm(node10, 2);
                 node11 = 3/5*node(elem(i, 2), :) + 2/5*node(elem(i, 1), :);
                 node11 = node11 /norm(node11, 2);
                 node12 = 4/5*node(elem(i, 2), :) + 1/5*node(elem(i, 1), :);
                 node12 = node12 /norm(node12, 2);
                 
                 node13 = 1/5*node(elem(i, 2), :) + 3/5*node(elem(i, 1), :) + 1/5*node(elem(i, 3), :);
                 node13 = node13 /norm(node13, 2);
                 
                 node14 = 2/5*node(elem(i, 2), :) + 2/5*node(elem(i, 1), :) + 1/5*node(elem(i, 3), :);
                 node14 = node14 /norm(node14, 2);
              

                 node15 = 3/5*node(elem(i, 2), :) + 1/5*node(elem(i, 1), :) + 1/5*node(elem(i, 3), :);
                 node15 = node15 /norm(node15, 2);

                 node16 = 1/5*node(elem(i, 2), :) + 2/5*node(elem(i, 1), :) + 2/5*node(elem(i, 3), :);
                 node16 = node16 /norm(node16, 2);

                 node17 = 2/5*node(elem(i, 2), :) + 1/5*node(elem(i, 1), :) + 2/5*node(elem(i, 3), :);
                 node17 = node17 /norm(node17, 2);

                 node18 = 1/5*node(elem(i, 2), :) + 1/5*node(elem(i, 1), :) + 3/5*node(elem(i, 3), :);
                 node18 = node18 /norm(node18, 2);


                 node_temp(4:21, :) = [node1;
                                     node2;
                                     node3;
                                     node4;
                                     node5;
                                     node6;
                                     node7;
                                     node8;
                                     node9;
                                     node10;
                                     node11;
                                     node12;
                                     node13;
                                     node14;
                                     node15;
                                     node16;
                                     node17;
                                     node18];




        for m=1:21
            %Derivatives of interpolation of geometry for building metric
            compx1 = compx1 + matrix(m, 1) * node_temp(m, 1);
            compy1 = compy1 + matrix(m, 1) * node_temp(m, 2);
            compz1 = compz1 + matrix(m, 1) * node_temp(m, 3);
            
            compx2 = compx2 + matrix(m, 2) * node_temp(m, 1);
            compy2 = compy2 + matrix(m, 2) * node_temp(m, 2);            
            compz2 = compz2 + matrix(m, 2) * node_temp(m, 3);

            hess_compx1 = hess_compx1 + hess_matrix(m, 1) * node_temp(m, 1);
            hess_compy1 = hess_compy1 + hess_matrix(m, 1) * node_temp(m, 2);
            hess_compz1 = hess_compz1 + hess_matrix(m, 1) * node_temp(m, 3);

            hess_compx2 = hess_compx2 + hess_matrix(m, 2) * node_temp(m, 1);
            hess_compy2 = hess_compy2 + hess_matrix(m, 2) * node_temp(m, 2);
            hess_compz2 = hess_compz2 + hess_matrix(m, 2) * node_temp(m, 3);

            hess_compx3 = hess_compx3 + hess_matrix(m, 3) * node_temp(m, 1);
            hess_compy3 = hess_compy3 + hess_matrix(m, 3) * node_temp(m, 2);
            hess_compz3 = hess_compz3 + hess_matrix(m, 3) * node_temp(m, 3);

            hess_compx4 = hess_compx4 + hess_matrix(m, 4) * node_temp(m, 1);
            hess_compy4 = hess_compy4 + hess_matrix(m, 4) * node_temp(m, 2);
            hess_compz4 = hess_compz4 + hess_matrix(m, 4) * node_temp(m, 3);

            point_comp1 = point_comp1 + point_matrix(m, 1) * node_temp(m, 1);
            point_comp2 = point_comp2 + point_matrix(m, 1) * node_temp(m, 2);
            point_comp3 = point_comp3 + point_matrix(m, 1) * node_temp(m, 3);
        end

        grad_ak_elem = [compx1, compx2; compy1 compy2; compz1 compz2]; 
        hess_ak_elem1 = [hess_compx1, hess_compx2; hess_compx3,hess_compx4];
        hess_ak_elem2 = [hess_compy1, hess_compy2; hess_compy3,hess_compy4];
        hess_ak_elem3 = [hess_compz1, hess_compz2; hess_compz3,hess_compz4];
        point_val = [point_comp1, point_comp2, point_comp3];



end        


end



