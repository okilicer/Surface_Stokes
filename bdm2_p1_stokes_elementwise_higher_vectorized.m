
clear 
%tic
%% Generate initial mesh

%  surfacedata = spheresurface_stokes();
% [node, elem] = surfacedata.initmesh();
% node = surfacedata.project(node);
% [node,elem] = optsurfacemesh(node,elem,surfacedata);
% [elem, idx, jac_det] = fixorder_surface(node, elem);
% NumRefinements = 1;
% if NumRefinements > 0
%     for i = 1:NumRefinements
%         [node, elem] = smeshuniformrefine(node,elem);
%         node = surfacedata.project(node);
%     end
% end
% [elem, idx, jac_det] = fixorder_surface(node, elem); 
% recquad = getquad_2d(6);
% a = 1.1;
% b = 1.2;
% c = 1.3;
a = 1;
b = 1;
c = 1;
tic
aa = [a b c];

node = [a, 0, 0; 0, b, 0; -a, 0, 0; 0, -b, 0; 0, 0, c; 0, 0, -c];

elem = [6, 1, 2; 6, 2, 3; 6, 3, 4; 6, 4, 1; 5, 1, 4; 5, 3, 4; 5, 3, 2; 5, 2, 1];

surf_deg = 3;

[elem, ~, ~] = fixorder_surface(node, elem);
surfacedata = my_ellipsoidsurface(a, b, c);
node = surfacedata.project(node);


NumRefinements = 6;

if NumRefinements > 0
    for i = 1:NumRefinements
        [node,elem] = smeshuniformrefine(node, elem);
        node = surfacedata.project(node);
    end
end
[elem, ~, ~] = fixorder_surface(node, elem); 

recquad=getquad_2d(4);


%%
% figure(1)
%     showmesh(node, elem);
% hold o
numel = size(elem, 1);

elemdof = zeros(numel, 12);
elemdofsign = zeros(numel, 12);
[elem2edge, edge, elem2edgeSign, edgeSign] = dofedge(elem);
T = auxstructure(elem);
edge2elem = T.edge2elem;

numedge = size(edge, 1);
    


%The following is for dof numbering and sign
degree = 3;
k = 1;
for i=1:numel
    while k < 4
        for j=1:3
            if elem2edgeSign(i, j) == -1
                elemdof(i, k) = elem2edge(i, j) * degree;
                elemdof(i, k+1) = elemdof(i, k) - 1;
                elemdof(i, k+2) = elemdof(i, k+1) - 1;
                %elemdof(i, k+3) = elemdof(i, k+2) - 1;
                elemdofsign(i, k) = elem2edgeSign(i, ceil(k/degree));
                elemdofsign(i, k+1) = elem2edgeSign(i, ceil((k+1)/degree));
                elemdofsign(i, k+2) = elem2edgeSign(i, ceil((k+2)/degree));
                %elemdofsign(i, k+3) = elem2edgeSign(i, ceil((k+3)/degree));
                %k = k + 3;
                k = k + degree;
            else
                elemdof(i, k) = elem2edge(i, j) * (degree) - degree + 1;
                elemdof(i, k+1) = elemdof(i, k) + 1;
                elemdof(i, k+2) = elemdof(i, k+1) + 1;
                %elemdof(i, k+3) = elemdof(i, k+2) + 1;
                elemdofsign(i, k) = elem2edgeSign(i, ceil(k/degree));
                elemdofsign(i, k+1) = elem2edgeSign(i, ceil((k+1)/degree));
                elemdofsign(i, k+2) = elem2edgeSign(i, ceil((k+2)/degree));
                %elemdofsign(i, k+3) = elem2edgeSign(i, ceil((k+3)/degree));
                %k = k + 3;
                k = k + degree;
            end
         end        
       
    end
    k = 1;
end
vector = (numedge*3+1):(numedge*3+numel*3);
elemdof(:, 10:12) = reshape(vector, numel, 3);
elemdofsign(:, 10:12) = 1;

numel = size(elem, 1);
numedge = size(edge, 1);


%Points

 noden = zeros(numel, 3, recquad.nqpts);
grad_ak = zeros(3, 2, recquad.nqpts, numel);

jac_det_elem = zeros(1, recquad.nqpts, numel);
nu_h = zeros(1, 3, recquad.nqpts, numel);
hessian_matrix = zeros(2, 2, recquad.nqpts, numel, 3);
summation = 0;
inverse = zeros(2, 3, recquad.nqpts, numel);
Ph_elem = zeros(3, 3, recquad.nqpts, numel);
for i=1:numel
    for m = 1:recquad.nqpts
        [point_val, h1, h2, h3, grad] = get_grad_ak(recquad.qpts(m, :), i, surf_deg, node, elem, surfacedata);

        grad_ak(:, :, m, i) = grad;
        jac_det_elem(:, m, i) = sqrt(det(grad_ak(:, :, m, i)'*grad_ak(:, :, m, i)));
        nu_h(:, :, m, i) = mycross(grad_ak(:, 1, m, i)', grad_ak(:, 2, m, i)');
        nu_h(:, :, m, i) = nu_h(:, :, m, i) / norm(nu_h(:, :, m, i), 2);
        pseudo_inverse_comp = grad_ak(:, :, m, i)' * grad_ak(:, :, m, i);
        inverse_comp = 1 / det(pseudo_inverse_comp)* [pseudo_inverse_comp(2, 2), ...
            -pseudo_inverse_comp(1, 2);
            -pseudo_inverse_comp(2, 1), pseudo_inverse_comp(1, 1)];
        hessian_matrix(:, :, m, i, 1) = h1;
        hessian_matrix(:, :, m, i, 2) = h2;
        hessian_matrix(:, :, m, i, 3) = h3;
        inverse(:, :, m, i) = inverse_comp * grad_ak(:, :, m, i)';
        noden(i, :, m) = point_val;
        Ph_elem(:, :, m, i) = eye(3) - nu_h(:, :, m, i)' * nu_h(:, :, m, i);
    end
end


h = max(max(sqrt(jac_det_elem/2)))

vector2 = 1:3*numel;
press_dof = reshape(vector2, numel, 3);
grhs = zeros(3*numel, 1);

iter = 0;
    for k=1:3
        kdofind = press_dof(:, k);
        num = zeros(numel, 1);
        for m = 1:recquad.nqpts
            p = noden(:, :, m);
            jac_det_vector = jac_det_elem(:, m, :);
            jac_det_vector = jac_det_vector(:);
            proj_point_surface = surfacedata.project(p);
            p1_basis = p1basis_orsan(recquad.qpts(m, :), k) ./ jac_det_vector;
            num = num + divusurf_el(proj_point_surface, aa).* p1_basis .* jac_det_vector * recquad.qwts(m);  
        end
        pkindex(iter*numel+1:(iter+1)*numel) = kdofind;
        grhs_vec(iter*numel+1:(iter+1)*numel) = num;
        iter = iter + 1;
    end

grhs = sparse(pkindex, 1, grhs_vec, 3*numel, 1);
% gavg = 2 * sum(grhs) / sum(sum(jac_det_elem));
% 
% grhs = grhs - gavg * jac_det / 2;


 
rvec = zeros(3*numedge+3*numel, 1);

iter = 0;
jindex = zeros(1, 1);

for j=1:12
    jdofind = elemdof(:, j);
    jdofsign = elemdofsign(:, j);
    rhs1 = zeros(numel, 1);
    
    for m=1:recquad.nqpts
        p = noden(:, :, m);
        proj_point_surface = surfacedata.project(p);
        nu = surfacedata.unitoutnormal(proj_point_surface);
        jac_det_vector = jac_det_elem(:, m, :);
        jac_det_vector = jac_det_vector(:);
        fval = -2 * divgradusurf_el(proj_point_surface, aa) + usurf_perp_el(proj_point_surface, aa) + gradpsurf_el(proj_point_surface, aa);
        nu_h_vec1 = nu_h(1, 1, m, :);
        nu_h_vec1 = nu_h_vec1(:);
        nu_h_vec2 = nu_h(1, 2, m, :);
        nu_h_vec2= nu_h_vec2(:);
        nu_h_vec3 = nu_h(1, 3, m, :);
        nu_h_vec3= nu_h_vec3(:);
        nu_h_vec = [nu_h_vec1, nu_h_vec2, nu_h_vec3];
        nuh_dot_fval = sum(fval.*nu_h_vec, 2);
        nu_dot_nuh = sum(nu.*nu_h_vec, 2);
        a1 = sum(grad_ak(1, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), j));
        a1 = jdofsign.*reshape(a1, [numel, 1])./jac_det_vector;
        b1 = sum(grad_ak(2, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), j));
        b1 = jdofsign.*reshape(b1, [numel, 1])./jac_det_vector;
        c1 = sum(grad_ak(3, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), j));
        c1 = jdofsign.*reshape(c1, [numel, 1])./jac_det_vector;
        vj1ref = [a1 b1 c1];
        fval = fval - nu .* nuh_dot_fval ./ nu_dot_nuh;
        rhs1 = rhs1+sum(fval.*vj1ref, 2).*jac_det_vector * recquad.qwts(m);
    end
    rvec(iter*numel+1:(iter+1)*numel) = rhs1;
    jindex(iter*numel+1:(iter+1)*numel) = jdofind;
    iter = iter + 1;
end

rhs = sparse(jindex, 1, rvec, 3*numedge+3*numel, 1);

jindex = zeros(1, 1);
kindex = zeros(1, 1);
Mvec = zeros(1, 1);
iter = 0;
%Calculation of the mass matrix


    for k=1:12
    kdofind = elemdof(:, k);
    kdofsign = elemdofsign(:, k);
    for j=1:12
        jdofind = elemdof(:, j);
        jdofsign = elemdofsign(:, j);
        
        Mjk = 0;
        for m = 1:recquad.nqpts
            jac_det_vector = jac_det_elem(:, m, :);
            jac_det_vector = jac_det_vector(:);
            a1 = sum(grad_ak(1, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), j));
            a1 = jdofsign.*reshape(a1, [numel, 1])./jac_det_vector;
            b1 = sum(grad_ak(2, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), j));
            b1 = jdofsign.*reshape(b1, [numel, 1])./jac_det_vector;
            c1 = sum(grad_ak(3, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), j));
            c1 = jdofsign.*reshape(c1, [numel, 1])./jac_det_vector;
            vj = [a1 b1 c1];
            a1 = sum(grad_ak(1, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), k));
            a1 = kdofsign.*reshape(a1, [numel, 1])./jac_det_vector;
            b1 = sum(grad_ak(2, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), k));
            b1 = kdofsign.*reshape(b1, [numel, 1])./jac_det_vector;
            c1 = sum(grad_ak(3, :, m, :) .* bdm_2basis_orsan(recquad.qpts(m, :), k));
            c1 = kdofsign.*reshape(c1, [numel, 1])./jac_det_vector;
            vk = [a1 b1 c1];

            Mjk = Mjk + sum(vk.*vj, 2) .*...
                 recquad.qwts(m) .* jac_det_vector; 
        end
        jindex(iter*numel+1:(iter+1)*numel) = jdofind;
        kindex(iter*numel+1:(iter+1)*numel) = kdofind;
        Mvec(iter*numel+1:(iter+1)*numel) = Mjk;
        iter = iter + 1;
    end
    end

M = sparse(jindex, kindex, Mvec, 3*numedge+3*numel, 3*numedge+3*numel);



Bjindex = zeros(3*numedge+3*numel, 1);
Bkindex = zeros(3*numel, 1);


iter = 0;

    for j = 1:12
        jdofind = elemdof(:, j);
        jdofsign = elemdofsign(:, j);
        for k=1:3
            kdofind = press_dof(:, k);
            div_term = zeros(numel, 1);
            for m = 1:recquad.nqpts
                jac_det_vector = jac_det_elem(:, m, :);
                jac_det_vector = jac_det_vector(:);
                first_term = p1basis_orsan(recquad.qpts(m, :), k)./jac_det_vector;
                trace_term = jdofsign .* trace(bdm2grad_basis_orsan(recquad.qpts(m, :), j))./jac_det_vector;
                div_term = div_term +  jac_det_vector ...
                       .* recquad.qwts(m) .* sum(first_term.*trace_term, 2);
          
            end
            Bjindex(iter*numel+1:(iter+1)*numel) = jdofind;
            Bkindex(iter*numel+1:(iter+1)*numel) = kdofind;
            Bvec(iter*numel+1:(iter+1)*numel) = div_term;
            iter = iter + 1;
        end

    end



B = sparse(Bjindex, Bkindex, Bvec, 3*numedge+3*numel, 3*numel);








Ajindex = zeros(144*numel, 1);
Akindex = zeros(144*numel, 1);
Avec = zeros(144*numel, 1);
iter = 0;
Hess_ak = zeros(6, 2, recquad.nqpts, numel);

    for j=1:12
        jdofind = elemdof(:, j);
        jdofsign = elemdofsign(:, j);
        for k=1:12
            kdofind = elemdof(:, k);
            kdofsign = elemdofsign(:, k);
            int11 = 0;
            for m=1:recquad.nqpts
                jac_det_vector = jac_det_elem(:, m, :);
                jac_det_vector = jac_det_vector(:);
                jacob1 = hessian_matrix(:, :, m, :, 1);
                jacob2 = hessian_matrix(:, :, m, :, 2);
                jacob3 = hessian_matrix(:, :, m, :, 3);
                a111 = jacob1(1, 1, :);
                a111 = reshape(a111, [numel, 1]);
                a112 = jacob1(1, 2, :);
                a112 = reshape(a112, [numel, 1]);
                a121 = a112;
                a122 = jacob1(2, 2, :);
                clear jacob1
                a122 = reshape(a122, [numel, 1]);
                a211 = jacob2(1, 1, :);
                a211 = reshape(a211, [numel, 1]);
                a212 = jacob2(1, 2, :);
                a212 = reshape(a212, [numel, 1]);
                a221 = a212;
                a222 = jacob2(2, 2, :);
                clear jacob2
                a222 = reshape(a222, [numel, 1]);
                a311 = jacob3(1, 1, :);
                a311 = reshape(a311, [numel, 1]);
                a312 = jacob3(1, 2, :);
                a312 = reshape(a312, [numel, 1]);
                a321 = a312;
                a322 = jacob3(2, 2, :);
                clear jacob3
                a322 = reshape(a322, [numel, 1]);
                a11 = grad_ak(1, 1, m, :);
                a11 = reshape(a11, [numel, 1]);
                a12 = grad_ak(1, 2, m, :);
                a12 = reshape(a12, [numel, 1]);
                a21 = grad_ak(2, 1, m, :);
                a21 = reshape(a21, [numel, 1]);
                a22 = grad_ak(2, 2, m, :);
                a22 = reshape(a22, [numel, 1]);
                a31 = grad_ak(3, 1, m, :);
                a31 = reshape(a31, [numel, 1]);
                a32 = grad_ak(3, 2, m, :);
                a32 = reshape(a32, [numel, 1]);
                jacob_matrixj1_11 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a111 a121], 2)./jac_det_vector;
                jacob_matrixj1_12 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a112 a122], 2)./jac_det_vector;
                jacob_matrixj1_21 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a211 a221], 2)./jac_det_vector;
                jacob_matrixj1_22 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a212 a222], 2)./jac_det_vector;
                jacob_matrixj1_31 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a311 a321], 2)./jac_det_vector;
                jacob_matrixj1_32 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a312 a322], 2)./jac_det_vector;
                jacob_matrix_j1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                   jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];


%                
                jacob_matrixj1_11 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a111 a121], 2)./jac_det_vector;
                jacob_matrixj1_12 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a112 a122], 2)./jac_det_vector;
                jacob_matrixj1_21 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a211 a221], 2)./jac_det_vector;
                jacob_matrixj1_22 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a212 a222], 2)./jac_det_vector;
                jacob_matrixj1_31 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a311 a321], 2)./jac_det_vector;
                jacob_matrixj1_32 =  sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a312 a322], 2)./jac_det_vector;
               
                jacob_matrix_k1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                     jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];


                grad_basis_termj = bdm2grad_basis_orsan(recquad.qpts(m, :), j);
                grad_basis_termk = bdm2grad_basis_orsan(recquad.qpts(m, :), k);
                
                
                
               
                jacob_matrixj1_11 = sum(grad_basis_termj(:, 1)' .* [a11 a12], 2)./jac_det_vector;
                jacob_matrixj1_12 = sum(grad_basis_termj(:, 2)' .* [a11 a12], 2)./jac_det_vector;
                jacob_matrixj1_21 = sum(grad_basis_termj(:, 1)' .* [a21 a22], 2)./jac_det_vector;
                jacob_matrixj1_22 = sum(grad_basis_termj(:, 2)' .* [a21 a22], 2)./jac_det_vector;
                jacob_matrixj1_31 = sum(grad_basis_termj(:, 1)' .* [a31 a32], 2)./jac_det_vector;
                jacob_matrixj1_32 = sum(grad_basis_termj(:, 2)' .* [a31 a32], 2)./jac_det_vector; 
                
                
                gradient_matrixj1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                       jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                
                sum_two_j1_1 = jacob_matrix_j1_1 + gradient_matrixj1_1;


                jacob_matrixj1_11 = sum(grad_basis_termk(:, 1)' .* [a11 a12], 2)./jac_det_vector;
                jacob_matrixj1_12 = sum(grad_basis_termk(:, 2)' .* [a11 a12], 2)./jac_det_vector;
                jacob_matrixj1_21 = sum(grad_basis_termk(:, 1)' .* [a21 a22], 2)./jac_det_vector;
                jacob_matrixj1_22 = sum(grad_basis_termk(:, 2)' .* [a21 a22], 2)./jac_det_vector;
                jacob_matrixj1_31 = sum(grad_basis_termk(:, 1)' .* [a31 a32], 2)./jac_det_vector;
                jacob_matrixj1_32 = sum(grad_basis_termk(:, 2)' .* [a31 a32], 2)./jac_det_vector; 
                
                gradient_matrixk1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                       jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                sum_two_k1_1 = jacob_matrix_k1_1 + gradient_matrixk1_1;

                det_grad_part1 = 2*a111.*a11.*a22.^2 + 2*a121.*a12.*a21.^2 + 2*a111.*a11.*a32.^2 + 2*a211.*a12.^2.*a21 + 2*a121.*a12.*a31.^2 + 2*a221.*a11.^2.*a22 + 2*a211.*a21.*a32.^2 + 2*a311.*a12.^2.*a31 + 2*a221.*a22.*a31.^2 + 2*a321.*a11.^2.*a32 + 2*a311.*a22.^2.*a31 + 2*a321.*a21.^2.*a32 - 2*a111.*a12.*a21.*a22 - 2*a121.*a11.*a21.*a22 - 2*a211.*a11.*a12.*a22 - 2*a221.*a11.*a12.*a21 - 2*a111.*a12.*a31.*a32 - 2*a121.*a11.*a31.*a32 - 2*a311.*a11.*a12.*a32 - 2*a321.*a11.*a12.*a31 - 2*a211.*a22.*a31.*a32 - 2*a221.*a21.*a31.*a32 - 2*a311.*a21.*a22.*a32 - 2*a321.*a21.*a22.*a31;
                det_grad_part2 = 2*a112.*a11.*a22.^2 + 2*a122.*a12.*a21.^2 + 2*a112.*a11.*a32.^2 + 2*a212.*a12.^2.*a21 + 2*a122.*a12.*a31.^2 + 2*a222.*a11.^2.*a22 + 2*a212.*a21.*a32.^2 + 2*a312.*a12.^2.*a31 + 2*a222.*a22.*a31.^2 + 2*a322.*a11.^2.*a32 + 2*a312.*a22.^2.*a31 + 2*a322.*a21.^2.*a32 - 2*a112.*a12.*a21.*a22 - 2*a122.*a11.*a21.*a22 - 2*a212.*a11.*a12.*a22 - 2*a222.*a11.*a12.*a21 - 2*a112.*a12.*a31.*a32 - 2*a122.*a11.*a31.*a32 - 2*a312.*a11.*a12.*a32 - 2*a322.*a11.*a12.*a31 - 2*a212.*a22.*a31.*a32 - 2*a222.*a21.*a31.*a32 - 2*a312.*a21.*a22.*a32 - 2*a322.*a21.*a22.*a31;
                
                detpartj1_11 = sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a11 a12], 2);
                detpartj1_21 = sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a21 a22], 2);
                detpartj1_31 = sum(bdm_2basis_orsan(recquad.qpts(m, :), j) .* [a31 a32], 2);
                
                detpartk1_11 = sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a11 a12], 2);
                detpartk1_21 = sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a21 a22], 2);
                detpartk1_31 = sum(bdm_2basis_orsan(recquad.qpts(m, :), k) .* [a31 a32], 2);
                
                jacob_matrixj1_11 = -1/2 * detpartj1_11 .* det_grad_part1./jac_det_vector.^3;
                jacob_matrixj1_12 = -1/2 * detpartj1_11 .* det_grad_part2./jac_det_vector.^3;
                jacob_matrixj1_21 = -1/2 * detpartj1_21 .* det_grad_part1./jac_det_vector.^3;
                jacob_matrixj1_22 = -1/2 * detpartj1_21 .* det_grad_part2./jac_det_vector.^3;
                jacob_matrixj1_31 = -1/2 * detpartj1_31 .* det_grad_part1./jac_det_vector.^3;
                jacob_matrixj1_32 = -1/2 * detpartj1_31 .* det_grad_part2./jac_det_vector.^3;
                
                det_grad_part_j1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                    jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                sum_three_j1 = sum_two_j1_1 + det_grad_part_j1;
                
                

                jacob_matrixj1_11 = -1/2 * detpartk1_11 .* det_grad_part1./jac_det_vector.^3;
                jacob_matrixj1_12 = -1/2 * detpartk1_11 .* det_grad_part2./jac_det_vector.^3;
                jacob_matrixj1_21 = -1/2 * detpartk1_21 .* det_grad_part1./jac_det_vector.^3;
                jacob_matrixj1_22 = -1/2 * detpartk1_21 .* det_grad_part2./jac_det_vector.^3;
                jacob_matrixj1_31 = -1/2 * detpartk1_31 .* det_grad_part1./jac_det_vector.^3;
                jacob_matrixj1_32 = -1/2 * detpartk1_31 .* det_grad_part2./jac_det_vector.^3;
                
                det_grad_part_k1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                    jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                sum_three_k1 = sum_two_k1_1 + det_grad_part_k1;
                

                 b11 = inverse(1, 1, m, :);
                b11= b11(:);
                b12 = inverse(1, 2, m, :);
                b12 = b12(:);
                b13 = inverse(1, 3, m, :);
                b13 = b13(:);
                b21 = inverse(2, 1, m, :);
                b21 = b21(:);
                b22 = inverse(2, 2, m, :);
                b22 = b22(:);
                b23 = inverse(2, 3, m, :);
                b23 = b23(:);

                grad_11 = sum_three_j1(:, 1).* b11 +  sum_three_j1(:, 2).* b21;  
                grad_12 = sum_three_j1(:, 1).* b12 +  sum_three_j1(:, 2).* b22;
                grad_13 = sum_three_j1(:, 1).* b13 +  sum_three_j1(:, 2).* b23; 
                grad_21 = sum_three_j1(:, 3).* b11 +  sum_three_j1(:, 4).* b21;  
                grad_22 = sum_three_j1(:, 3).* b12 +  sum_three_j1(:, 4).* b22;
                grad_23 = sum_three_j1(:, 3).* b13 +  sum_three_j1(:, 4).* b23;
                grad_31 = sum_three_j1(:, 5).* b11 +  sum_three_j1(:, 6).* b21;  
                grad_32 = sum_three_j1(:, 5).* b12 +  sum_three_j1(:, 6).* b22;
                grad_33 = sum_three_j1(:, 5).* b13 +  sum_three_j1(:, 6).* b23;
                
                 clear sum_three_j1

                grad_j1 = [grad_11, grad_12, grad_13, grad_21, grad_22, ...
                           grad_23, grad_31, grad_32, grad_33];

                grad_11 = sum_three_k1(:, 1).* b11 +  sum_three_k1(:, 2).* b21;  
                grad_12 = sum_three_k1(:, 1).* b12 +  sum_three_k1(:, 2).* b22;
                grad_13 = sum_three_k1(:, 1).* b13 +  sum_three_k1(:, 2).* b23; 
                grad_21 = sum_three_k1(:, 3).* b11 +  sum_three_k1(:, 4).* b21;  
                grad_22 = sum_three_k1(:, 3).* b12 +  sum_three_k1(:, 4).* b22;
                grad_23 = sum_three_k1(:, 3).* b13 +  sum_three_k1(:, 4).* b23;
                grad_31 = sum_three_k1(:, 5).* b11 +  sum_three_k1(:, 6).* b21;  
                grad_32 = sum_three_k1(:, 5).* b12 +  sum_three_k1(:, 6).* b22;
                grad_33 = sum_three_k1(:, 5).* b13 +  sum_three_k1(:, 6).* b23;
                
                clear sum_three_k1

                grad_k1 = [grad_11, grad_12, grad_13, grad_21, grad_22, ...
                           grad_23, grad_31, grad_32, grad_33];


                def_j1 = [grad_j1(:, 1), 1/2*(grad_j1(:, 2)+grad_j1(:, 4)),...
                           1/2*(grad_j1(:, 3)+grad_j1(:, 7)), 1/2*(grad_j1(:, 2)+grad_j1(:, 4)),...
                           grad_j1(:, 5), 1/2*(grad_j1(:, 6)+grad_j1(:, 8)), ...
                           1/2*(grad_j1(:, 3)+grad_j1(:, 7)), 1/2*(grad_j1(:, 6)+grad_j1(:, 8)),...
                           grad_j1(:, 9)];
                 clear grad_j1

                
                def_k1 = [grad_k1(:, 1), 1/2*(grad_k1(:, 2)+grad_k1(:, 4)),...
                           1/2*(grad_k1(:, 3)+grad_k1(:, 7)), 1/2*(grad_k1(:, 2)+grad_k1(:, 4)),...
                           grad_k1(:, 5), 1/2*(grad_k1(:, 6)+grad_k1(:, 8)), ...
                           1/2*(grad_k1(:, 3)+grad_k1(:, 7)), 1/2*(grad_k1(:, 6)+grad_k1(:, 8)),...
                           grad_k1(:, 9)];
                 clear grad_k1






                 b11 = Ph_elem(1, 1, m, :);
                b11 = b11(:);
                b12 = Ph_elem(1, 2, m, :);
                b12 = b12(:);
                b13 = Ph_elem(1, 3, m, :);
                b13 = b13(:);
                b21 = Ph_elem(2, 1, m, :);
                b21 = b21(:);
                b22 = Ph_elem(2, 2, m, :);
                b22 = b22(:);
                b23 = Ph_elem(2, 3, m, :);
                b23 = b23(:);
                b31 = Ph_elem(3, 1, m, :);
                b31 = b31(:);
                b32 = Ph_elem(3, 2, m, :);
                b32 = b32(:);
                b33 = Ph_elem(3, 3, m, :);
                b33 = b33(:);

                ph_left_11 = (def_j1(:, 1) .* b11) + (def_j1(:, 4) .* b12) + (def_j1(:, 7) .* b13);
                ph_left_12 = (def_j1(:, 2) .* b11) + (def_j1(:, 5) .* b12) + (def_j1(:, 8) .* b13);
                ph_left_13 = (def_j1(:, 3) .* b11) + (def_j1(:, 6) .* b12) + (def_j1(:, 9) .* b13);
                ph_left_21 = (def_j1(:, 1) .* b21) + (def_j1(:, 4) .* b22) + (def_j1(:, 7) .* b23);
                ph_left_22 = (def_j1(:, 2) .* b21) + (def_j1(:, 5) .* b22) + (def_j1(:, 8) .* b23);
                ph_left_23 = (def_j1(:, 3) .* b21) + (def_j1(:, 6) .* b22) + (def_j1(:, 9) .* b23);
                ph_left_31 = (def_j1(:, 1) .* b31) + (def_j1(:, 4) .* b32) + (def_j1(:, 7) .* b33);
                ph_left_32 = (def_j1(:, 2) .* b31) + (def_j1(:, 5) .* b32) + (def_j1(:, 8) .* b33);
                ph_left_33 = (def_j1(:, 3) .* b31) + (def_j1(:, 6) .* b32) + (def_j1(:, 9) .* b33);
 
                ph_left_def_j1 = [ph_left_11, ph_left_12, ph_left_13,...
                                  ph_left_21, ph_left_22, ph_left_23,...
                                  ph_left_31, ph_left_32, ph_left_33];
                clear def_j1





                ph_left_11 = (def_k1(:, 1) .* b11) + (def_k1(:, 4) .* b12) + (def_k1(:, 7) .* b13);
                ph_left_12 = (def_k1(:, 2) .* b11) + (def_k1(:, 5) .* b12) + (def_k1(:, 8) .* b13);
                ph_left_13 = (def_k1(:, 3) .* b11) + (def_k1(:, 6) .* b12) + (def_k1(:, 9) .* b13);
                ph_left_21 = (def_k1(:, 1) .* b21) + (def_k1(:, 4) .* b22) + (def_k1(:, 7) .* b23);
                ph_left_22 = (def_k1(:, 2) .* b21) + (def_k1(:, 5) .* b22) + (def_k1(:, 8) .* b23);
                ph_left_23 = (def_k1(:, 3) .* b21) + (def_k1(:, 6) .* b22) + (def_k1(:, 9) .* b23);
                ph_left_31 = (def_k1(:, 1) .* b31) + (def_k1(:, 4) .* b32) + (def_k1(:, 7) .* b33);
                ph_left_32 = (def_k1(:, 2) .* b31) + (def_k1(:, 5) .* b32) + (def_k1(:, 8) .* b33);
                ph_left_33 = (def_k1(:, 3) .* b31) + (def_k1(:, 6) .* b32) + (def_k1(:, 9) .* b33);
 
                ph_left_def_k1 = [ph_left_11, ph_left_12, ph_left_13,...
                                  ph_left_21, ph_left_22, ph_left_23,...
                                  ph_left_31, ph_left_32, ph_left_33];
                 clear def_k1

                ph_right_11 = (ph_left_def_j1(:, 1) .* b11) + (ph_left_def_j1(:, 2) .* b21) + (ph_left_def_j1(:, 3) .* b31);
                ph_right_12 = (ph_left_def_j1(:, 1) .* b12) + (ph_left_def_j1(:, 2) .* b22) + (ph_left_def_j1(:, 3) .* b32);
                ph_right_13 = (ph_left_def_j1(:, 1) .* b13) + (ph_left_def_j1(:, 2) .* b23) + (ph_left_def_j1(:, 3) .* b33);
                ph_right_21 = (ph_left_def_j1(:, 4) .* b11) + (ph_left_def_j1(:, 5) .* b21) + (ph_left_def_j1(:, 6) .* b31);
                ph_right_22 = (ph_left_def_j1(:, 4) .* b12) + (ph_left_def_j1(:, 5) .* b22) + (ph_left_def_j1(:, 6) .* b32);
                ph_right_23 = (ph_left_def_j1(:, 4) .* b13) + (ph_left_def_j1(:, 5) .* b23) + (ph_left_def_j1(:, 6) .* b33);
                ph_right_31 = (ph_left_def_j1(:, 7) .* b11) + (ph_left_def_j1(:, 8) .* b21) + (ph_left_def_j1(:, 9) .* b31);
                ph_right_32 = (ph_left_def_j1(:, 7) .* b12) + (ph_left_def_j1(:, 8) .* b22) + (ph_left_def_j1(:, 9) .* b32);
                ph_right_33 = (ph_left_def_j1(:, 7) .* b13) + (ph_left_def_j1(:, 8) .* b23) + (ph_left_def_j1(:, 9) .* b33);
 
                full_def_j_1 = [ph_right_11, ph_right_12, ph_right_13,...
                               ph_right_21, ph_right_22, ph_right_23,...
                               ph_right_31, ph_right_32, ph_right_33];
                
                 clear ph_left_def_j1




                ph_right_11 = (ph_left_def_k1(:, 1) .* b11) + (ph_left_def_k1(:, 2) .* b21) + (ph_left_def_k1(:, 3) .* b31);
                ph_right_12 = (ph_left_def_k1(:, 1) .* b12) + (ph_left_def_k1(:, 2) .* b22) + (ph_left_def_k1(:, 3) .* b32);
                ph_right_13 = (ph_left_def_k1(:, 1) .* b13) + (ph_left_def_k1(:, 2) .* b23) + (ph_left_def_k1(:, 3) .* b33);
                ph_right_21 = (ph_left_def_k1(:, 4) .* b11) + (ph_left_def_k1(:, 5) .* b21) + (ph_left_def_k1(:, 6) .* b31);
                ph_right_22 = (ph_left_def_k1(:, 4) .* b12) + (ph_left_def_k1(:, 5) .* b22) + (ph_left_def_k1(:, 6) .* b32);
                ph_right_23 = (ph_left_def_k1(:, 4) .* b13) + (ph_left_def_k1(:, 5) .* b23) + (ph_left_def_k1(:, 6) .* b33);
                ph_right_31 = (ph_left_def_k1(:, 7) .* b11) + (ph_left_def_k1(:, 8) .* b21) + (ph_left_def_k1(:, 9) .* b31);
                ph_right_32 = (ph_left_def_k1(:, 7) .* b12) + (ph_left_def_k1(:, 8) .* b22) + (ph_left_def_k1(:, 9) .* b32);
                ph_right_33 = (ph_left_def_k1(:, 7) .* b13) + (ph_left_def_k1(:, 8) .* b23) + (ph_left_def_k1(:, 9) .* b33);
 
                full_def_k_1 = [ph_right_11, ph_right_12, ph_right_13,...
                               ph_right_21, ph_right_22, ph_right_23,...
                               ph_right_31, ph_right_32, ph_right_33];
                
                 clear ph_left_def_k1


                int11 = int11 + jdofsign.*kdofsign.*sum(full_def_j_1 .* full_def_k_1, 2) * recquad.qwts(m).*jac_det_vector;

            end
            %Compute stiffness matrix (volume terms)
            Ajindex(iter*numel+1:(iter+1)*numel) = jdofind;
            Akindex(iter*numel+1:(iter+1)*numel) = kdofind;
            Avec(iter*numel+1:(iter+1)*numel) = int11;
            iter = iter + 1;
        end
    end  


% 

toc
edgequad = getquad_1d(4);

beta = 10;
% 
AAjindex = zeros(numedge*edgequad.nqpts*12, 1);
AAkindex = zeros(numedge*edgequad.nqpts*12, 1);
AAvec = zeros(numedge*edgequad.nqpts*12, 1);




grad_ak_edge = zeros(3, 2, edgequad.nqpts, numedge, 2);
jac_det_elem_edge = zeros(1, edgequad.nqpts, numedge, 2);
nu_h_edge = zeros(1, 3, edgequad.nqpts, numedge, 2);
hessian_matrix_edge = zeros(2, 2, edgequad.nqpts, numedge, 3, 2);
inverse_edge = zeros(2, 3, edgequad.nqpts, numedge, 2);
Ph_elem_edge = zeros(3, 3, edgequad.nqpts, numedge, 2);
edge_tangent = zeros(1, 3, edgequad.nqpts, numedge);
h_edge = zeros(1, edgequad.nqpts, numedge);
t_conormal_edge = zeros(1, 3, edgequad.nqpts, numedge, 2);
vector_basis_edge = zeros(1, 2, edgequad.nqpts, numedge, 12, 2);
grad_basis_edge = zeros(2, 2, edgequad.nqpts, numedge, 12, 2);
for i=1:numedge
    edgesigni = edgeSign(i);
    localedge1 = edge2elem(i, 3);
    localedge2 = edge2elem(i, 4);
    elem1 = edge2elem(i, 1);
    elem2 = edge2elem(i, 2);
    for m = 1:edgequad.nqpts
        
        lambda1 = zeros(1, 3);
        lambda2 = lambda1;
        if elem2edgeSign(elem1, localedge1)==1
            lambda1(1, mod(localedge1, 3)+1) = edgequad.qpts(m, 1);
            lambda1(1, mod(localedge1+1, 3)+1) = edgequad.qpts(m, 2);
        else
            lambda1(1, mod(localedge1, 3)+1) = edgequad.qpts(m, 2);
            lambda1(1, mod(localedge1+1, 3)+1) = edgequad.qpts(m, 1);
        end
        if elem2edgeSign(elem2, localedge2)==1
            lambda2(1, mod(localedge2, 3)+1) = edgequad.qpts(m, 1);
            lambda2(1, mod(localedge2+1, 3)+1) = edgequad.qpts(m, 2);
        else
            lambda2(1, mod(localedge2, 3)+1) = edgequad.qpts(m, 2);
            lambda2(1, mod(localedge2+1, 3)+1) = edgequad.qpts(m, 1);
        end
        [~, jacob11, jacob21, jacob31, grad_ak1] = get_grad_ak(lambda1, elem1, surf_deg, node, elem, surfacedata);
        [~, jacob12, jacob22, jacob32, grad_ak2] = get_grad_ak(lambda2, elem2, surf_deg, node, elem, surfacedata);

        grad_ak_edge(:, :, m, i, 1) = grad_ak1;
        grad_ak_edge(:, :, m, i, 2) = grad_ak2;
        jac_det_elem_edge(:, m, i, 1) = sqrt(det(grad_ak1'*grad_ak1));
        jac_det_elem_edge(:, m, i, 2) = sqrt(det(grad_ak2'*grad_ak2));
        nu_h_edge(:, :, m, i, 1) = mycross(grad_ak_edge(:, 1, m, i, 1)', grad_ak_edge(:, 2, m, i, 1)');
        nu_h_edge(:, :, m, i, 1) = nu_h_edge(:, :, m, i, 1) / norm(nu_h_edge(:, :, m, i, 1), 2);
        nu_h_edge(:, :, m, i, 2) = mycross(grad_ak_edge(:, 1, m, i, 2)', grad_ak_edge(:, 2, m, i, 2)');
        nu_h_edge(:, :, m, i, 2) = nu_h_edge(:, :, m, i, 2) / norm(nu_h_edge(:, :, m, i, 2), 2);
        
        pseudo_inverse_comp1 = grad_ak1' * grad_ak1;
        pseudo_inverse_comp2 = grad_ak2' * grad_ak2;
        inverse_comp1 = 1 / det(pseudo_inverse_comp1)* [pseudo_inverse_comp1(2, 2), ...
            -pseudo_inverse_comp1(1, 2);
            -pseudo_inverse_comp1(2, 1), pseudo_inverse_comp1(1, 1)];
        inverse_comp2 = 1 / det(pseudo_inverse_comp2)* [pseudo_inverse_comp2(2, 2), ...
            -pseudo_inverse_comp2(1, 2);
            -pseudo_inverse_comp2(2, 1), pseudo_inverse_comp2(1, 1)];
        hessian_matrix_edge(:, :, m, i, 1, 1) = jacob11;
        hessian_matrix_edge(:, :, m, i, 2, 1) = jacob21;
        hessian_matrix_edge(:, :, m, i, 3, 1) = jacob31;
        hessian_matrix_edge(:, :, m, i, 1, 2) = jacob12;
        hessian_matrix_edge(:, :, m, i, 2, 2) = jacob22;
        hessian_matrix_edge(:, :, m, i, 3, 2) = jacob32;
        inverse_edge(:, :, m, i, 1) = inverse_comp1 * grad_ak1';
        inverse_edge(:, :, m, i, 2) = inverse_comp2 * grad_ak2';
        Ph_elem_edge(:, :, m, i, 1) = eye(3) - nu_h_edge(:, :, m, i, 1)' * nu_h_edge(:, :, m, i, 1);
        Ph_elem_edge(:, :, m, i, 2) = eye(3) - nu_h_edge(:, :, m, i, 2)' * nu_h_edge(:, :, m, i, 2);
        if (localedge1 == 1) && (elem2edgeSign(elem1, localedge1)==-1)
            edge_tangent(:, :, m, i) = grad_ak1 * [-1; 1];
        elseif (localedge1 == 1) && (elem2edgeSign(elem1, localedge1)==1)
            edge_tangent(:, :, m, i) = -1 * grad_ak1 * [-1; 1];
        elseif (localedge1 == 2) && (elem2edgeSign(elem1, localedge1)==1)
            edge_tangent(:, :, m, i) = grad_ak1 * [0; 1];
        elseif (localedge1 == 2) && (elem2edgeSign(elem1, localedge1)==-1)
            edge_tangent(:, :, m, i) = -1 * grad_ak1 * [0; 1];
        elseif(localedge1 == 3) && (elem2edgeSign(elem1, localedge1)==1)
            edge_tangent(:, :, m, i) = -1 * grad_ak1 * [1; 0];
        elseif(localedge1 == 3) && (elem2edgeSign(elem1, localedge1)==-1)
            edge_tangent(:, :, m, i) = grad_ak1 * [1; 0];
        end

        h_edge(:, m, i) = sqrt(sum(edge_tangent(:, :, m, i)*edge_tangent(:, :, m, i)', 2));
        edge_tangent(:, :, m, i) = edge_tangent(:, :, m, i) / norm(edge_tangent(:, :, m, i), 2);
        if (localedge1 == localedge2) && (edgesigni==-1)
            t_conormal_edge(:, :, m, i, 1) = -mycross(nu_h_edge(:, :, m, i, 1), edge_tangent(:, :, m, i));
            t_conormal_edge(:, :, m, i, 2) = mycross(nu_h_edge(:, :, m, i, 2), edge_tangent(:, :, m, i));
            
        elseif(localedge1 ~= localedge2) && (edgesigni==-1)
            t_conormal_edge(:, :, m, i, 1) = -mycross(nu_h_edge(:, :, m, i, 1), edge_tangent(:, :, m, i));
            t_conormal_edge(:, :, m, i, 2) = mycross(nu_h_edge(:, :, m, i, 2), edge_tangent(:, :, m, i));
        else
            t_conormal_edge(:, :, m, i, 1) = mycross(nu_h_edge(:, :, m, i, 1), edge_tangent(:, :, m, i));
            t_conormal_edge(:, :, m, i, 2) = -mycross(nu_h_edge(:, :, m, i, 2), edge_tangent(:, :, m, i));
        end
        for j=1:12
            vector_basis_edge(:, :, m, i, j, 1) = bdm_2basis_orsan(lambda1, j);
            vector_basis_edge(:, :, m, i, j, 2) = bdm_2basis_orsan(lambda2, j);
            grad_basis_edge(:, :, m, i, j, 1) = bdm2grad_basis_orsan(lambda1, j);
            grad_basis_edge(:, :, m, i, j, 2) = bdm2grad_basis_orsan(lambda2, j);
        end
    end
end





tic
iter = 0;
    for j=1:12
        jds1 = elemdofsign(edge2elem(:, 1), j);
        jds2 = elemdofsign(edge2elem(:, 2), j);

        for k=1:12
            kds1 = elemdofsign(edge2elem(:, 1), k);
            kds2 = elemdofsign(edge2elem(:, 2), k);
            A1jk = 0;  
            A2jk = 0;
            A3jk = 0;
            A4jk = 0;
            for m=1:edgequad.nqpts
                jac_det_vector1 = jac_det_elem_edge(:, m, :, 1);
                jac_det_vector1 = jac_det_vector1(:);
                jac_det_vector2 = jac_det_elem_edge(:, m, :, 2);
                jac_det_vector2 = jac_det_vector2(:);

                jacob11 = hessian_matrix_edge(:, :, m, :, 1, 1);
                jacob21 = hessian_matrix_edge(:, :, m, :, 2, 1);
                jacob31 = hessian_matrix_edge(:, :, m, :, 3, 1);
                a111 = jacob11(1, 1, :);
                a111 = reshape(a111, [numedge, 1]);
                a112 = jacob11(1, 2, :);
                a112 = reshape(a112, [numedge, 1]);
                a121 = a112;
                a122 = jacob11(2, 2, :);
                a122 = reshape(a122, [numedge, 1]);
                a211 = jacob21(1, 1, :);
                a211 = reshape(a211, [numedge, 1]);
                a212 = jacob21(1, 2, :);
                a212 = reshape(a212, [numedge, 1]);
                a221 = a212;
                a222 = jacob21(2, 2, :);
                a222 = reshape(a222, [numedge, 1]);
                a311 = jacob31(1, 1, :);
                a311 = reshape(a311, [numedge, 1]);
                a312 = jacob31(1, 2, :);
                a312 = reshape(a312, [numedge, 1]);
                a321 = a312;
                a322 = jacob31(2, 2, :);
                a322 = reshape(a322, [numedge, 1]);
                a11 = grad_ak_edge(1, 1, m, :, 1);
                a11 = reshape(a11, [numedge, 1]);
                a12 = grad_ak_edge(1, 2, m, :, 1);
                a12 = reshape(a12, [numedge, 1]);
                a21 = grad_ak_edge(2, 1, m, :, 1);
                a21 = reshape(a21, [numedge, 1]);
                a22 = grad_ak_edge(2, 2, m, :, 1);
                a22 = reshape(a22, [numedge, 1]);
                a31 = grad_ak_edge(3, 1, m, :, 1);
                a31 = reshape(a31, [numedge, 1]);
                a32 = grad_ak_edge(3, 2, m, :, 1);
                a32 = reshape(a32, [numedge, 1]);

                jacob12 = hessian_matrix_edge(:, :, m, :, 1, 2);
                jacob22 = hessian_matrix_edge(:, :, m, :, 2, 2);
                jacob32 = hessian_matrix_edge(:, :, m, :, 3, 2);
                b111 = jacob12(1, 1, :);
                b111 = reshape(b111, [numedge, 1]);
                b112 = jacob12(1, 2, :);
                b112 = reshape(b112, [numedge, 1]);
                b121 = b112;
                b122 = jacob12(2, 2, :);
                b122 = reshape(b122, [numedge, 1]);
                b211 = jacob22(1, 1, :);
                b211 = reshape(b211, [numedge, 1]);
                b212 = jacob22(1, 2, :);
                b212 = reshape(b212, [numedge, 1]);
                b221 = b212;
                b222 = jacob22(2, 2, :);
                b222 = reshape(b222, [numedge, 1]);
                b311 = jacob32(1, 1, :);
                b311 = reshape(b311, [numedge, 1]);
                b312 = jacob32(1, 2, :);
                b312 = reshape(b312, [numedge, 1]);
                b321 = b312;
                b322 = jacob32(2, 2, :);
                b322 = reshape(b322, [numedge, 1]);

                b11 = grad_ak_edge(1, 1, m, :, 2);
                b11 = reshape(b11, [numedge, 1]);
                b12 = grad_ak_edge(1, 2, m, :, 2);
                b12 = reshape(b12, [numedge, 1]);
                b21 = grad_ak_edge(2, 1, m, :, 2);
                b21 = reshape(b21, [numedge, 1]);
                b22 = grad_ak_edge(2, 2, m, :, 2);
                b22 = reshape(b22, [numedge, 1]);
                b31 = grad_ak_edge(3, 1, m, :, 2);
                b31 = reshape(b31, [numedge, 1]);
                b32 = grad_ak_edge(3, 2, m, :, 2);
                b32 = reshape(b32, [numedge, 1]);




      
                
                





                jac_det_elem1 = jac_det_elem_edge(:, m, :, 1);
                jac_det_elem1 = jac_det_elem1(:);
                jac_det_elem2 = jac_det_elem_edge(:, m, :, 2);
                jac_det_elem2 = jac_det_elem2(:);

                basis11j = vector_basis_edge(1, 1, m, :, j, 1);
                basis11j = reshape(basis11j, [numedge, 1]);
                basis12j = vector_basis_edge(1, 2, m, :, j, 1);
                basis12j = reshape(basis12j, [numedge, 1]);
            
                basis21j = vector_basis_edge(1, 1, m, :, j, 2);
                basis21j = reshape(basis21j, [numedge, 1]);
                basis22j = vector_basis_edge(1, 2, m, :, j, 2);
                basis22j = reshape(basis22j, [numedge, 1]);

                basis11k = vector_basis_edge(1, 1, m, :, k, 1);
                basis11k = reshape(basis11k, [numedge, 1]);
                basis12k = vector_basis_edge(1, 2, m, :, k, 1);
                basis12k = reshape(basis12k, [numedge, 1]);
            
                basis21k = vector_basis_edge(1, 1, m, :, k, 2);
                basis21k = reshape(basis21k, [numedge, 1]);
                basis22k = vector_basis_edge(1, 2, m, :, k, 2);
                basis22k = reshape(basis22k, [numedge, 1]);



                edgetangent_1 = edge_tangent(1, 1, m, :);
                edgetangent_1 = reshape(edgetangent_1, [numedge, 1]);
                edgetangent_2 = edge_tangent(1, 2, m, :);
                edgetangent_2 = reshape(edgetangent_2, [numedge, 1]);
                edgetangent_3 = edge_tangent(1, 3, m, :);
                edgetangent_3 = reshape(edgetangent_3, [numedge, 1]);
                
                uplus1 = sum([a11 a12] .* [basis11j basis12j], 2);
                uplus2 = sum([a21 a22] .* [basis11j basis12j], 2);
                uplus3 = sum([a31 a32] .* [basis11j basis12j], 2);

                u1j = sum([uplus1 uplus2 uplus3] .* [edgetangent_1 edgetangent_2 edgetangent_3], 2)./jac_det_elem1;
                
%                 u1j = edgetangent1 * uplus;


                uminus1 = sum([b11 b12] .* [basis21j basis22j], 2);
                uminus2 = sum([b21 b22] .* [basis21j basis22j], 2);
                uminus3 = sum([b31 b32] .* [basis21j basis22j], 2);
                
                u2j = sum([uminus1 uminus2 uminus3] .* [edgetangent_1 edgetangent_2 edgetangent_3], 2)./jac_det_elem2;
%                 u2j = edgetangent1 * uminus; 

                uplus1 = sum([a11 a12] .* [basis11k basis12k], 2);
                uplus2 = sum([a21 a22] .* [basis11k basis12k], 2);
                uplus3 = sum([a31 a32] .* [basis11k basis12k], 2);

%                 u1k = edgetangent1 * uplus;
                u1k = sum([uplus1 uplus2 uplus3] .* [edgetangent_1 edgetangent_2 edgetangent_3], 2)./jac_det_elem1;
               
                uminus1 = sum([b11 b12] .* [basis21k basis22k], 2);
                uminus2 = sum([b21 b22] .* [basis21k basis22k], 2);
                uminus3 = sum([b31 b32] .* [basis21k basis22k], 2);
                
                u2k = sum([uminus1 uminus2 uminus3] .* [edgetangent_1 edgetangent_2 edgetangent_3], 2)./jac_det_elem2;
%                 u2k = edgetangent1 * uminus; 
                
                
                jacob_matrixj1_11 =  sum([basis11j basis12j] .* [a111 a121], 2)./jac_det_elem1;
                jacob_matrixj1_12 =  sum([basis11j basis12j] .* [a112 a122], 2)./jac_det_elem1;
                jacob_matrixj1_21 =  sum([basis11j basis12j] .* [a211 a221], 2)./jac_det_elem1;
                jacob_matrixj1_22 =  sum([basis11j basis12j] .* [a212 a222], 2)./jac_det_elem1;
                jacob_matrixj1_31 =  sum([basis11j basis12j] .* [a311 a321], 2)./jac_det_elem1;
                jacob_matrixj1_32 =  sum([basis11j basis12j] .* [a312 a322], 2)./jac_det_elem1;
                jacob_matrix_j1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                   jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];

                jacob_matrixj1_11 =  sum([basis11k basis12k] .* [a111 a121], 2)./jac_det_elem1;
                jacob_matrixj1_12 =  sum([basis11k basis12k] .* [a112 a122], 2)./jac_det_elem1;
                jacob_matrixj1_21 =  sum([basis11k basis12k] .* [a211 a221], 2)./jac_det_elem1;
                jacob_matrixj1_22 =  sum([basis11k basis12k] .* [a212 a222], 2)./jac_det_elem1;
                jacob_matrixj1_31 =  sum([basis11k basis12k] .* [a311 a321], 2)./jac_det_elem1;
                jacob_matrixj1_32 =  sum([basis11k basis12k] .* [a312 a322], 2)./jac_det_elem1;
                jacob_matrix_k1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                   jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];


                jacob_matrixj1_11 =  sum([basis21j basis22j] .* [b111 b121], 2)./jac_det_elem2;
                jacob_matrixj1_12 =  sum([basis21j basis22j] .* [b112 b122], 2)./jac_det_elem2;
                jacob_matrixj1_21 =  sum([basis21j basis22j] .* [b211 b221], 2)./jac_det_elem2;
                jacob_matrixj1_22 =  sum([basis21j basis22j] .* [b212 b222], 2)./jac_det_elem2;
                jacob_matrixj1_31 =  sum([basis21j basis22j] .* [b311 b321], 2)./jac_det_elem2;
                jacob_matrixj1_32 =  sum([basis21j basis22j] .* [b312 b322], 2)./jac_det_elem2;
                jacob_matrix_j1_2 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                   jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];


                jacob_matrixj1_11 =  sum([basis21k basis22k] .* [b111 b121], 2)./jac_det_elem2;
                jacob_matrixj1_12 =  sum([basis21k basis22k] .* [b112 b122], 2)./jac_det_elem2;
                jacob_matrixj1_21 =  sum([basis21k basis22k] .* [b211 b221], 2)./jac_det_elem2;
                jacob_matrixj1_22 =  sum([basis21k basis22k] .* [b212 b222], 2)./jac_det_elem2;
                jacob_matrixj1_31 =  sum([basis21k basis22k] .* [b311 b321], 2)./jac_det_elem2;
                jacob_matrixj1_32 =  sum([basis21k basis22k] .* [b312 b322], 2)./jac_det_elem2;
                jacob_matrix_k1_2 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                   jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
 

                grad_basistermj1 = grad_basis_edge(:, :, m, :, j, 1);
                grad_basistermj2 = grad_basis_edge(:, :, m, :, j, 2);
                grad_basistermk1 = grad_basis_edge(:, :, m, :, k, 1);
                grad_basistermk2 = grad_basis_edge(:, :, m, :, k, 2);
                grad_basistermj11 = [reshape(grad_basistermj1(1:1, 1:1, :, :), [numedge, 1]) reshape(grad_basistermj1(2:2, 1:1, :, :), [numedge, 1])] ;
                grad_basistermj12 = [reshape(grad_basistermj1(1:1, 2:2, :, :), [numedge, 1]) reshape(grad_basistermj1(2:2, 2:2, :, :), [numedge, 1])] ;
                grad_basistermj21 = [reshape(grad_basistermj2(1:1, 1:1, :, :), [numedge, 1]) reshape(grad_basistermj2(2:2, 1:1, :, :), [numedge, 1])] ;
                grad_basistermj22 = [reshape(grad_basistermj2(1:1, 2:2, :, :), [numedge, 1]) reshape(grad_basistermj2(2:2, 2:2, :, :), [numedge, 1])] ;
                grad_basistermk11 = [reshape(grad_basistermk1(1:1, 1:1, :, :), [numedge, 1]) reshape(grad_basistermk1(2:2, 1:1, :, :), [numedge, 1])] ;
                grad_basistermk12 = [reshape(grad_basistermk1(1:1, 2:2, :, :), [numedge, 1]) reshape(grad_basistermk1(2:2, 2:2, :, :), [numedge, 1])] ;
                grad_basistermk21 = [reshape(grad_basistermk2(1:1, 1:1, :, :), [numedge, 1]) reshape(grad_basistermk2(2:2, 1:1, :, :), [numedge, 1])] ;
                grad_basistermk22 = [reshape(grad_basistermk2(1:1, 2:2, :, :), [numedge, 1]) reshape(grad_basistermk2(2:2, 2:2, :, :), [numedge, 1])] ;




                jacob_matrixj1_11 = sum(grad_basistermj11 .* [a11 a12], 2)./jac_det_elem1;
                jacob_matrixj1_12 = sum(grad_basistermj12 .* [a11 a12], 2)./jac_det_elem1;
                jacob_matrixj1_21 = sum(grad_basistermj11 .* [a21 a22], 2)./jac_det_elem1;
                jacob_matrixj1_22 = sum(grad_basistermj12 .* [a21 a22], 2)./jac_det_elem1;
                jacob_matrixj1_31 = sum(grad_basistermj11 .* [a31 a32], 2)./jac_det_elem1;
                jacob_matrixj1_32 = sum(grad_basistermj12 .* [a31 a32], 2)./jac_det_elem1; 
                
                
                gradient_matrixj1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                       jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];

%                 grad_basis_edge(:, :, m, i, j, 1)
                sum_two_j1_1 = jacob_matrix_j1_1 + gradient_matrixj1_1;


                jacob_matrixj1_11 = sum(grad_basistermj21 .* [b11 b12], 2)./jac_det_elem2;
                jacob_matrixj1_12 = sum(grad_basistermj22 .* [b11 b12], 2)./jac_det_elem2;
                jacob_matrixj1_21 = sum(grad_basistermj21 .* [b21 b22], 2)./jac_det_elem2;
                jacob_matrixj1_22 = sum(grad_basistermj22 .* [b21 b22], 2)./jac_det_elem2;
                jacob_matrixj1_31 = sum(grad_basistermj21 .* [b31 b32], 2)./jac_det_elem2;
                jacob_matrixj1_32 = sum(grad_basistermj22 .* [b31 b32], 2)./jac_det_elem2; 
                
                
                gradient_matrixj1_2 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                       jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];

%                 grad_basis_edge(:, :, m, i, j, 1)
                sum_two_j1_2 = jacob_matrix_j1_2 + gradient_matrixj1_2;


                jacob_matrixj1_11 = sum(grad_basistermk11 .* [a11 a12], 2)./jac_det_elem1;
                jacob_matrixj1_12 = sum(grad_basistermk12 .* [a11 a12], 2)./jac_det_elem1;
                jacob_matrixj1_21 = sum(grad_basistermk11 .* [a21 a22], 2)./jac_det_elem1;
                jacob_matrixj1_22 = sum(grad_basistermk12 .* [a21 a22], 2)./jac_det_elem1;
                jacob_matrixj1_31 = sum(grad_basistermk11 .* [a31 a32], 2)./jac_det_elem1;
                jacob_matrixj1_32 = sum(grad_basistermk12 .* [a31 a32], 2)./jac_det_elem1; 
                
                gradient_matrixk1_1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                       jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                sum_two_k1_1 = jacob_matrix_k1_1 + gradient_matrixk1_1;




                jacob_matrixj1_11 = sum(grad_basistermk21 .* [b11 b12], 2)./jac_det_elem2;
                jacob_matrixj1_12 = sum(grad_basistermk22 .* [b11 b12], 2)./jac_det_elem2;
                jacob_matrixj1_21 = sum(grad_basistermk21 .* [b21 b22], 2)./jac_det_elem2;
                jacob_matrixj1_22 = sum(grad_basistermk22 .* [b21 b22], 2)./jac_det_elem2;
                jacob_matrixj1_31 = sum(grad_basistermk21 .* [b31 b32], 2)./jac_det_elem2;
                jacob_matrixj1_32 = sum(grad_basistermk22 .* [b31 b32], 2)./jac_det_elem2; 
                
                
                gradient_matrixk1_2 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                       jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];

%                 grad_basis_edge(:, :, m, i, j, 1)
                sum_two_k1_2 = jacob_matrix_k1_2 + gradient_matrixk1_2;
                
                
                
                
                
                
                
                
                
                


                
                
                
                
                
                
                
                %det_grad_part11an = jacob11(1, 1)*grad_ak1(2, 2)^2*grad_ak1(1, 1) + grad_ak1(2, 2)^2*grad_ak1(1, 1)*jacob11(1, 1) + jacob11(2, 1)*grad_ak1(2, 1)^2*grad_ak1(1, 2) + grad_ak1(2, 1)^2*grad_ak1(1, 2)*jacob11(2, 1) + jacob11(1, 1)*grad_ak1(3, 2)^2*grad_ak1(1, 1) + jacob21(1, 1)*grad_ak1(1, 2)^2*grad_ak1(2, 1) + grad_ak1(1, 2)^2*grad_ak1(2, 1)*jacob21(1, 1) + grad_ak1(3, 2)^2*grad_ak1(1, 1)*jacob11(1, 1) + jacob11(2, 1)*grad_ak1(3, 1)^2*grad_ak1(1, 2) + jacob21(2, 1)*grad_ak1(1, 1)^2*grad_ak1(2, 2) + grad_ak1(1, 1)^2*grad_ak1(2, 2)*jacob21(2, 1) + grad_ak1(3, 1)^2*grad_ak1(1, 2)*jacob11(2, 1) + jacob21(1, 1)*grad_ak1(3, 2)^2*grad_ak1(2, 1) + jacob31(1, 1)*grad_ak1(1, 2)^2*grad_ak1(3, 1) + grad_ak1(1, 2)^2*grad_ak1(3, 1)*jacob31(1, 1) + grad_ak1(3, 2)^2*grad_ak1(2, 1)*jacob21(1, 1) + jacob21(2, 1)*grad_ak1(3, 1)^2*grad_ak1(2, 2) + jacob31(2, 1)*grad_ak1(1, 1)^2*grad_ak1(3, 2) + grad_ak1(1, 1)^2*grad_ak1(3, 2)*jacob31(2, 1) + grad_ak1(3, 1)^2*grad_ak1(2, 2)*jacob21(2, 1) + jacob31(1, 1)*grad_ak1(2, 2)^2*grad_ak1(3, 1) + grad_ak1(2, 2)^2*grad_ak1(3, 1)*jacob31(1, 1) + jacob31(2, 1)*grad_ak1(2, 1)^2*grad_ak1(3, 2) + grad_ak1(2, 1)^2*grad_ak1(3, 2)*jacob31(2, 1) - jacob11(1, 1)*grad_ak1(2, 2)*grad_ak1(1, 2)*grad_ak1(2, 1) - jacob11(2, 1)*grad_ak1(2, 1)*grad_ak1(1, 1)*grad_ak1(2, 2) - jacob21(1, 1)*grad_ak1(1, 2)*grad_ak1(1, 1)*grad_ak1(2, 2) - jacob21(2, 1)*grad_ak1(1, 1)*grad_ak1(1, 2)*grad_ak1(2, 1) - grad_ak1(1, 1)*grad_ak1(2, 2)*grad_ak1(2, 1)*jacob11(2, 1) - grad_ak1(1, 1)*grad_ak1(2, 2)*grad_ak1(1, 2)*jacob21(1, 1) - grad_ak1(1, 2)*grad_ak1(2, 1)*grad_ak1(2, 2)*jacob11(1, 1) - grad_ak1(1, 2)*grad_ak1(2, 1)*grad_ak1(1, 1)*jacob21(2, 1) - jacob11(1, 1)*grad_ak1(3, 2)*grad_ak1(1, 2)*grad_ak1(3, 1) - jacob11(2, 1)*grad_ak1(3, 1)*grad_ak1(1, 1)*grad_ak1(3, 2) - jacob31(1, 1)*grad_ak1(1, 2)*grad_ak1(1, 1)*grad_ak1(3, 2) - jacob31(2, 1)*grad_ak1(1, 1)*grad_ak1(1, 2)*grad_ak1(3, 1) - grad_ak1(1, 1)*grad_ak1(3, 2)*grad_ak1(3, 1)*jacob11(2, 1) - grad_ak1(1, 1)*grad_ak1(3, 2)*grad_ak1(1, 2)*jacob31(1, 1) - grad_ak1(1, 2)*grad_ak1(3, 1)*grad_ak1(3, 2)*jacob11(1, 1) - grad_ak1(1, 2)*grad_ak1(3, 1)*grad_ak1(1, 1)*jacob31(2, 1) - jacob21(1, 1)*grad_ak1(3, 2)*grad_ak1(2, 2)*grad_ak1(3, 1) - jacob21(2, 1)*grad_ak1(3, 1)*grad_ak1(2, 1)*grad_ak1(3, 2) - jacob31(1, 1)*grad_ak1(2, 2)*grad_ak1(2, 1)*grad_ak1(3, 2) - jacob31(2, 1)*grad_ak1(2, 1)*grad_ak1(2, 2)*grad_ak1(3, 1) - grad_ak1(2, 1)*grad_ak1(3, 2)*grad_ak1(3, 1)*jacob21(2, 1) - grad_ak1(2, 1)*grad_ak1(3, 2)*grad_ak1(2, 2)*jacob31(1, 1) - grad_ak1(2, 2)*grad_ak1(3, 1)*grad_ak1(3, 2)*jacob21(1, 1) - grad_ak1(2, 2)*grad_ak1(3, 1)*grad_ak1(2, 1)*jacob31(2, 1);
                det_grad_part11 = 2*a111.*a11.*a22.^2 + 2*a121.*a12.*a21.^2 + 2*a111.*a11.*a32.^2 + 2*a211.*a12.^2.*a21 + 2*a121.*a12.*a31.^2 + 2*a221.*a11.^2.*a22 + 2*a211.*a21.*a32.^2 + 2*a311.*a12.^2.*a31 + 2*a221.*a22.*a31.^2 + 2*a321.*a11.^2.*a32 + 2*a311.*a22.^2.*a31 + 2*a321.*a21.^2.*a32 - 2*a111.*a12.*a21.*a22 - 2*a121.*a11.*a21.*a22 - 2*a211.*a11.*a12.*a22 - 2*a221.*a11.*a12.*a21 - 2*a111.*a12.*a31.*a32 - 2*a121.*a11.*a31.*a32 - 2*a311.*a11.*a12.*a32 - 2*a321.*a11.*a12.*a31 - 2*a211.*a22.*a31.*a32 - 2*a221.*a21.*a31.*a32 - 2*a311.*a21.*a22.*a32 - 2*a321.*a21.*a22.*a31;
                det_grad_part12 = 2*a112.*a11.*a22.^2 + 2*a122.*a12.*a21.^2 + 2*a112.*a11.*a32.^2 + 2*a212.*a12.^2.*a21 + 2*a122.*a12.*a31.^2 + 2*a222.*a11.^2.*a22 + 2*a212.*a21.*a32.^2 + 2*a312.*a12.^2.*a31 + 2*a222.*a22.*a31.^2 + 2*a322.*a11.^2.*a32 + 2*a312.*a22.^2.*a31 + 2*a322.*a21.^2.*a32 - 2*a112.*a12.*a21.*a22 - 2*a122.*a11.*a21.*a22 - 2*a212.*a11.*a12.*a22 - 2*a222.*a11.*a12.*a21 - 2*a112.*a12.*a31.*a32 - 2*a122.*a11.*a31.*a32 - 2*a312.*a11.*a12.*a32 - 2*a322.*a11.*a12.*a31 - 2*a212.*a22.*a31.*a32 - 2*a222.*a21.*a31.*a32 - 2*a312.*a21.*a22.*a32 - 2*a322.*a21.*a22.*a31;
                %det_grad_part12an = jacob11(1, 2)*grad_ak1(2, 2)^2*grad_ak1(1, 1) + grad_ak1(2, 2)^2*grad_ak1(1, 1)*jacob11(1, 2) + jacob11(2, 2)*grad_ak1(2, 1)^2*grad_ak1(1, 2) + grad_ak1(2, 1)^2*grad_ak1(1, 2)*jacob11(2, 2) + jacob11(1, 2)*grad_ak1(3, 2)^2*grad_ak1(1, 1) + jacob21(1, 2)*grad_ak1(1, 2)^2*grad_ak1(2, 1) + grad_ak1(1, 2)^2*grad_ak1(2, 1)*jacob21(1, 2) + grad_ak1(3, 2)^2*grad_ak1(1, 1)*jacob11(1, 2) + jacob11(2, 2)*grad_ak1(3, 1)^2*grad_ak1(1, 2) + jacob21(2, 2)*grad_ak1(1, 1)^2*grad_ak1(2, 2) + grad_ak1(1, 1)^2*grad_ak1(2, 2)*jacob21(2, 2) + grad_ak1(3, 1)^2*grad_ak1(1, 2)*jacob11(2, 2) + jacob21(1, 2)*grad_ak1(3, 2)^2*grad_ak1(2, 1) + jacob31(1, 2)*grad_ak1(1, 2)^2*grad_ak1(3, 1) + grad_ak1(1, 2)^2*grad_ak1(3, 1)*jacob31(1, 2) + grad_ak1(3, 2)^2*grad_ak1(2, 1)*jacob21(1, 2) + jacob21(2, 2)*grad_ak1(3, 1)^2*grad_ak1(2, 2) + jacob31(2, 2)*grad_ak1(1, 1)^2*grad_ak1(3, 2) + grad_ak1(1, 1)^2*grad_ak1(3, 2)*jacob31(2, 2) + grad_ak1(3, 1)^2*grad_ak1(2, 2)*jacob21(2, 2) + jacob31(1, 2)*grad_ak1(2, 2)^2*grad_ak1(3, 1) + grad_ak1(2, 2)^2*grad_ak1(3, 1)*jacob31(1, 2) + jacob31(2, 2)*grad_ak1(2, 1)^2*grad_ak1(3, 2) + grad_ak1(2, 1)^2*grad_ak1(3, 2)*jacob31(2, 2) - jacob11(1, 2)*grad_ak1(2, 2)*grad_ak1(1, 2)*grad_ak1(2, 1) - jacob11(2, 2)*grad_ak1(2, 1)*grad_ak1(1, 1)*grad_ak1(2, 2) - jacob21(1, 2)*grad_ak1(1, 2)*grad_ak1(1, 1)*grad_ak1(2, 2) - jacob21(2, 2)*grad_ak1(1, 1)*grad_ak1(1, 2)*grad_ak1(2, 1) - grad_ak1(1, 1)*grad_ak1(2, 2)*grad_ak1(2, 1)*jacob11(2, 2) - grad_ak1(1, 1)*grad_ak1(2, 2)*grad_ak1(1, 2)*jacob21(1, 2) - grad_ak1(1, 2)*grad_ak1(2, 1)*grad_ak1(2, 2)*jacob11(1, 2) - grad_ak1(1, 2)*grad_ak1(2, 1)*grad_ak1(1, 1)*jacob21(2, 2) - jacob11(1, 2)*grad_ak1(3, 2)*grad_ak1(1, 2)*grad_ak1(3, 1) - jacob11(2, 2)*grad_ak1(3, 1)*grad_ak1(1, 1)*grad_ak1(3, 2) - jacob31(1, 2)*grad_ak1(1, 2)*grad_ak1(1, 1)*grad_ak1(3, 2) - jacob31(2, 2)*grad_ak1(1, 1)*grad_ak1(1, 2)*grad_ak1(3, 1) - grad_ak1(1, 1)*grad_ak1(3, 2)*grad_ak1(3, 1)*jacob11(2, 2) - grad_ak1(1, 1)*grad_ak1(3, 2)*grad_ak1(1, 2)*jacob31(1, 2) - grad_ak1(1, 2)*grad_ak1(3, 1)*grad_ak1(3, 2)*jacob11(1, 2) - grad_ak1(1, 2)*grad_ak1(3, 1)*grad_ak1(1, 1)*jacob31(2, 2) - jacob21(1, 2)*grad_ak1(3, 2)*grad_ak1(2, 2)*grad_ak1(3, 1) - jacob21(2, 2)*grad_ak1(3, 1)*grad_ak1(2, 1)*grad_ak1(3, 2) - jacob31(1, 2)*grad_ak1(2, 2)*grad_ak1(2, 1)*grad_ak1(3, 2) - jacob31(2, 2)*grad_ak1(2, 1)*grad_ak1(2, 2)*grad_ak1(3, 1) - grad_ak1(2, 1)*grad_ak1(3, 2)*grad_ak1(3, 1)*jacob21(2, 2) - grad_ak1(2, 1)*grad_ak1(3, 2)*grad_ak1(2, 2)*jacob31(1, 2) - grad_ak1(2, 2)*grad_ak1(3, 1)*grad_ak1(3, 2)*jacob21(1, 2) - grad_ak1(2, 2)*grad_ak1(3, 1)*grad_ak1(2, 1)*jacob31(2, 2);

%                 
                det_grad11 = [det_grad_part11 det_grad_part12];


                det_grad_part21 = 2*b111.*b11.*b22.^2 + 2*b121.*b12.*b21.^2 + 2*b111.*b11.*b32.^2 + 2*b211.*b12.^2.*b21 + 2*b121.*b12.*b31.^2 + 2*b221.*b11.^2.*b22 + 2*b211.*b21.*b32.^2 + 2*b311.*b12.^2.*b31 + 2*b221.*b22.*b31.^2 + 2*b321.*b11.^2.*b32 + 2*b311.*b22.^2.*b31 + 2*b321.*b21.^2.*b32 - 2*b111.*b12.*b21.*b22 - 2*b121.*b11.*b21.*b22 - 2*b211.*b11.*b12.*b22 - 2*b221.*b11.*b12.*b21 - 2*b111.*b12.*b31.*b32 - 2*b121.*b11.*b31.*b32 - 2*b311.*b11.*b12.*b32 - 2*b321.*b11.*b12.*b31 - 2*b211.*b22.*b31.*b32 - 2*b221.*b21.*b31.*b32 - 2*b311.*b21.*b22.*b32 - 2*b321.*b21.*b22.*b31;
                det_grad_part22 = 2*b112.*b11.*b22.^2 + 2*b122.*b12.*b21.^2 + 2*b112.*b11.*b32.^2 + 2*b212.*b12.^2.*b21 + 2*b122.*b12.*b31.^2 + 2*b222.*b11.^2.*b22 + 2*b212.*b21.*b32.^2 + 2*b312.*b12.^2.*b31 + 2*b222.*b22.*b31.^2 + 2*b322.*b11.^2.*b32 + 2*b312.*b22.^2.*b31 + 2*b322.*b21.^2.*b32 - 2*b112.*b12.*b21.*b22 - 2*b122.*b11.*b21.*b22 - 2*b212.*b11.*b12.*b22 - 2*b222.*b11.*b12.*b21 - 2*b112.*b12.*b31.*b32 - 2*b122.*b11.*b31.*b32 - 2*b312.*b11.*b12.*b32 - 2*b322.*b11.*b12.*b31 - 2*b212.*b22.*b31.*b32 - 2*b222.*b21.*b31.*b32 - 2*b312.*b21.*b22.*b32 - 2*b322.*b21.*b22.*b31;
                det_grad12 = [det_grad_part21 det_grad_part22];



                detpartj1_11 = sum([basis11j basis12j] .* [a11 a12], 2);
                detpartj1_21 = sum([basis11j basis12j] .* [a21 a22], 2);
                detpartj1_31 = sum([basis11j basis12j] .* [a31 a32], 2);

                detpartj2_11 = sum([basis21j basis22j] .* [b11 b12], 2);
                detpartj2_21 = sum([basis21j basis22j] .* [b21 b22], 2);
                detpartj2_31 = sum([basis21j basis22j] .* [b31 b32], 2);
                
                detpartk1_11 = sum([basis11k basis12k] .* [a11 a12], 2);
                detpartk1_21 = sum([basis11k basis12k] .* [a21 a22], 2);
                detpartk1_31 = sum([basis11k basis12k] .* [a31 a32], 2);

                detpartk2_11 = sum([basis21k basis22k] .* [b11 b12], 2);
                detpartk2_21 = sum([basis21k basis22k] .* [b21 b22], 2);
                detpartk2_31 = sum([basis21k basis22k] .* [b31 b32], 2);
                
                jacob_matrixj1_11 = -1/2 * detpartj1_11 .* det_grad_part11./jac_det_elem1.^3;
                jacob_matrixj1_12 = -1/2 * detpartj1_11 .* det_grad_part12./jac_det_elem1.^3;
                jacob_matrixj1_21 = -1/2 * detpartj1_21 .* det_grad_part11./jac_det_elem1.^3;
                jacob_matrixj1_22 = -1/2 * detpartj1_21 .* det_grad_part12./jac_det_elem1.^3;
                jacob_matrixj1_31 = -1/2 * detpartj1_31 .* det_grad_part11./jac_det_elem1.^3;
                jacob_matrixj1_32 = -1/2 * detpartj1_31 .* det_grad_part12./jac_det_elem1.^3;
                
                det_grad_part_j1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                    jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                sum_three_j1 = sum_two_j1_1 + det_grad_part_j1;



                jacob_matrixj1_11 = -1/2 * detpartj2_11 .* det_grad_part21./jac_det_elem2.^3;
                jacob_matrixj1_12 = -1/2 * detpartj2_11 .* det_grad_part22./jac_det_elem2.^3;
                jacob_matrixj1_21 = -1/2 * detpartj2_21 .* det_grad_part21./jac_det_elem2.^3;
                jacob_matrixj1_22 = -1/2 * detpartj2_21 .* det_grad_part22./jac_det_elem2.^3;
                jacob_matrixj1_31 = -1/2 * detpartj2_31 .* det_grad_part21./jac_det_elem2.^3;
                jacob_matrixj1_32 = -1/2 * detpartj2_31 .* det_grad_part22./jac_det_elem2.^3;
                
                det_grad_part_j2 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                    jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                sum_three_j2 = sum_two_j1_2 + det_grad_part_j2;


                jacob_matrixj1_11 = -1/2 * detpartk1_11 .* det_grad_part11./jac_det_elem1.^3;
                jacob_matrixj1_12 = -1/2 * detpartk1_11 .* det_grad_part12./jac_det_elem1.^3;
                jacob_matrixj1_21 = -1/2 * detpartk1_21 .* det_grad_part11./jac_det_elem1.^3;
                jacob_matrixj1_22 = -1/2 * detpartk1_21 .* det_grad_part12./jac_det_elem1.^3;
                jacob_matrixj1_31 = -1/2 * detpartk1_31 .* det_grad_part11./jac_det_elem1.^3;
                jacob_matrixj1_32 = -1/2 * detpartk1_31 .* det_grad_part12./jac_det_elem1.^3;
                
                det_grad_part_k1 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                    jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                sum_three_k1 = sum_two_k1_1 + det_grad_part_k1;
 
                


                jacob_matrixj1_11 = -1/2 * detpartk2_11 .* det_grad_part21./jac_det_elem2.^3;
                jacob_matrixj1_12 = -1/2 * detpartk2_11 .* det_grad_part22./jac_det_elem2.^3;
                jacob_matrixj1_21 = -1/2 * detpartk2_21 .* det_grad_part21./jac_det_elem2.^3;
                jacob_matrixj1_22 = -1/2 * detpartk2_21 .* det_grad_part22./jac_det_elem2.^3;
                jacob_matrixj1_31 = -1/2 * detpartk2_31 .* det_grad_part21./jac_det_elem2.^3;
                jacob_matrixj1_32 = -1/2 * detpartk2_31 .* det_grad_part22./jac_det_elem2.^3;
                
                det_grad_part_k2 = [jacob_matrixj1_11, jacob_matrixj1_12, jacob_matrixj1_21, ...
                                    jacob_matrixj1_22, jacob_matrixj1_31, jacob_matrixj1_32];
                
                sum_three_k2 = sum_two_k1_2 + det_grad_part_k2;

                b11 = inverse_edge(1, 1, m, :, 1);
                b11= b11(:);
                b12 = inverse_edge(1, 2, m, :, 1);
                b12 = b12(:);
                b13 = inverse_edge(1, 3, m, :, 1);
                b13 = b13(:);
                b21 = inverse_edge(2, 1, m, :, 1);
                b21 = b21(:);
                b22 = inverse_edge(2, 2, m, :, 1);
                b22 = b22(:);
                b23 = inverse_edge(2, 3, m, :, 1);
                b23 = b23(:);


                grad_11 = sum_three_j1(:, 1).* b11 +  sum_three_j1(:, 2).* b21;  
                grad_12 = sum_three_j1(:, 1).* b12 +  sum_three_j1(:, 2).* b22;
                grad_13 = sum_three_j1(:, 1).* b13 +  sum_three_j1(:, 2).* b23; 
                grad_21 = sum_three_j1(:, 3).* b11 +  sum_three_j1(:, 4).* b21;  
                grad_22 = sum_three_j1(:, 3).* b12 +  sum_three_j1(:, 4).* b22;
                grad_23 = sum_three_j1(:, 3).* b13 +  sum_three_j1(:, 4).* b23;
                grad_31 = sum_three_j1(:, 5).* b11 +  sum_three_j1(:, 6).* b21;  
                grad_32 = sum_three_j1(:, 5).* b12 +  sum_three_j1(:, 6).* b22;
                grad_33 = sum_three_j1(:, 5).* b13 +  sum_three_j1(:, 6).* b23;

                grad_j1 = [grad_11, grad_12, grad_13, grad_21, grad_22, ...
                           grad_23, grad_31, grad_32, grad_33];


                grad_11 = sum_three_k1(:, 1).* b11 +  sum_three_k1(:, 2).* b21;  
                grad_12 = sum_three_k1(:, 1).* b12 +  sum_three_k1(:, 2).* b22;
                grad_13 = sum_three_k1(:, 1).* b13 +  sum_three_k1(:, 2).* b23; 
                grad_21 = sum_three_k1(:, 3).* b11 +  sum_three_k1(:, 4).* b21;  
                grad_22 = sum_three_k1(:, 3).* b12 +  sum_three_k1(:, 4).* b22;
                grad_23 = sum_three_k1(:, 3).* b13 +  sum_three_k1(:, 4).* b23;
                grad_31 = sum_three_k1(:, 5).* b11 +  sum_three_k1(:, 6).* b21;  
                grad_32 = sum_three_k1(:, 5).* b12 +  sum_three_k1(:, 6).* b22;
                grad_33 = sum_three_k1(:, 5).* b13 +  sum_three_k1(:, 6).* b23;

                grad_k1 = [grad_11, grad_12, grad_13, grad_21, grad_22, ...
                           grad_23, grad_31, grad_32, grad_33];



                b11 = inverse_edge(1, 1, m, :, 2);
                b11 = b11(:);
                b12 = inverse_edge(1, 2, m, :, 2);
                b12 = b12(:);
                b13 = inverse_edge(1, 3, m, :, 2);
                b13 = b13(:);
                b21 = inverse_edge(2, 1, m, :, 2);
                b21 = b21(:);
                b22 = inverse_edge(2, 2, m, :, 2);
                b22 = b22(:);
                b23 = inverse_edge(2, 3, m, :, 2);
                b23 = b23(:);



                grad_11 = sum_three_j2(:, 1).* b11 +  sum_three_j2(:, 2).* b21;  
                grad_12 = sum_three_j2(:, 1).* b12 +  sum_three_j2(:, 2).* b22;
                grad_13 = sum_three_j2(:, 1).* b13 +  sum_three_j2(:, 2).* b23; 
                grad_21 = sum_three_j2(:, 3).* b11 +  sum_three_j2(:, 4).* b21;  
                grad_22 = sum_three_j2(:, 3).* b12 +  sum_three_j2(:, 4).* b22;
                grad_23 = sum_three_j2(:, 3).* b13 +  sum_three_j2(:, 4).* b23;
                grad_31 = sum_three_j2(:, 5).* b11 +  sum_three_j2(:, 6).* b21;  
                grad_32 = sum_three_j2(:, 5).* b12 +  sum_three_j2(:, 6).* b22;
                grad_33 = sum_three_j2(:, 5).* b13 +  sum_three_j2(:, 6).* b23;

                grad_j2 = [grad_11, grad_12, grad_13, grad_21, grad_22, ...
                           grad_23, grad_31, grad_32, grad_33];


                grad_11 = sum_three_k2(:, 1).* b11 +  sum_three_k2(:, 2).* b21;  
                grad_12 = sum_three_k2(:, 1).* b12 +  sum_three_k2(:, 2).* b22;
                grad_13 = sum_three_k2(:, 1).* b13 +  sum_three_k2(:, 2).* b23; 
                grad_21 = sum_three_k2(:, 3).* b11 +  sum_three_k2(:, 4).* b21;  
                grad_22 = sum_three_k2(:, 3).* b12 +  sum_three_k2(:, 4).* b22;
                grad_23 = sum_three_k2(:, 3).* b13 +  sum_three_k2(:, 4).* b23;
                grad_31 = sum_three_k2(:, 5).* b11 +  sum_three_k2(:, 6).* b21;  
                grad_32 = sum_three_k2(:, 5).* b12 +  sum_three_k2(:, 6).* b22;
                grad_33 = sum_three_k2(:, 5).* b13 +  sum_three_k2(:, 6).* b23;

                grad_k2 = [grad_11, grad_12, grad_13, grad_21, grad_22, ...
                           grad_23, grad_31, grad_32, grad_33];

                def_j1 = [grad_j1(:, 1), 1/2*(grad_j1(:, 2)+grad_j1(:, 4)),...
                           1/2*(grad_j1(:, 3)+grad_j1(:, 7)), 1/2*(grad_j1(:, 2)+grad_j1(:, 4)),...
                           grad_j1(:, 5), 1/2*(grad_j1(:, 6)+grad_j1(:, 8)), ...
                           1/2*(grad_j1(:, 3)+grad_j1(:, 7)), 1/2*(grad_j1(:, 6)+grad_j1(:, 8)),...
                           grad_j1(:, 9)];


                def_k1 = [grad_k1(:, 1), 1/2*(grad_k1(:, 2)+grad_k1(:, 4)),...
                           1/2*(grad_k1(:, 3)+grad_k1(:, 7)), 1/2*(grad_k1(:, 2)+grad_k1(:, 4)),...
                           grad_k1(:, 5), 1/2*(grad_k1(:, 6)+grad_k1(:, 8)), ...
                           1/2*(grad_k1(:, 3)+grad_k1(:, 7)), 1/2*(grad_k1(:, 6)+grad_k1(:, 8)),...
                           grad_k1(:, 9)];


                def_j2 = [grad_j2(:, 1), 1/2*(grad_j2(:, 2)+grad_j2(:, 4)),...
                           1/2*(grad_j2(:, 3)+grad_j2(:, 7)), 1/2*(grad_j2(:, 2)+grad_j2(:, 4)),...
                           grad_j2(:, 5), 1/2*(grad_j2(:, 6)+grad_j2(:, 8)), ...
                           1/2*(grad_j2(:, 3)+grad_j2(:, 7)), 1/2*(grad_j2(:, 6)+grad_j2(:, 8)),...
                           grad_j2(:, 9)];


                def_k2 = [grad_k2(:, 1), 1/2*(grad_k2(:, 2)+grad_k2(:, 4)),...
                           1/2*(grad_k2(:, 3)+grad_k2(:, 7)), 1/2*(grad_k2(:, 2)+grad_k2(:, 4)),...
                           grad_k2(:, 5), 1/2*(grad_k2(:, 6)+grad_k2(:, 8)), ...
                           1/2*(grad_k2(:, 3)+grad_k2(:, 7)), 1/2*(grad_k2(:, 6)+grad_k2(:, 8)),...
                           grad_k2(:, 9)];





                b11 = Ph_elem_edge(1, 1, m, :, 1);
                b11 = b11(:);
                b12 = Ph_elem_edge(1, 2, m, :, 1);
                b12 = b12(:);
                b13 = Ph_elem_edge(1, 3, m, :, 1);
                b13 = b13(:);
                b21 = Ph_elem_edge(2, 1, m, :, 1);
                b21 = b21(:);
                b22 = Ph_elem_edge(2, 2, m, :, 1);
                b22 = b22(:);
                b23 = Ph_elem_edge(2, 3, m, :, 1);
                b23 = b23(:);
                b31 = Ph_elem_edge(3, 1, m, :, 1);
                b31 = b31(:);
                b32 = Ph_elem_edge(3, 2, m, :, 1);
                b32 = b32(:);
                b33 = Ph_elem_edge(3, 3, m, :, 1);
                b33 = b33(:);

                ph_left_11 = (def_j1(:, 1) .* b11) + (def_j1(:, 4) .* b12) + (def_j1(:, 7) .* b13);
                ph_left_12 = (def_j1(:, 2) .* b11) + (def_j1(:, 5) .* b12) + (def_j1(:, 8) .* b13);
                ph_left_13 = (def_j1(:, 3) .* b11) + (def_j1(:, 6) .* b12) + (def_j1(:, 9) .* b13);
                ph_left_21 = (def_j1(:, 1) .* b21) + (def_j1(:, 4) .* b22) + (def_j1(:, 7) .* b23);
                ph_left_22 = (def_j1(:, 2) .* b21) + (def_j1(:, 5) .* b22) + (def_j1(:, 8) .* b23);
                ph_left_23 = (def_j1(:, 3) .* b21) + (def_j1(:, 6) .* b22) + (def_j1(:, 9) .* b23);
                ph_left_31 = (def_j1(:, 1) .* b31) + (def_j1(:, 4) .* b32) + (def_j1(:, 7) .* b33);
                ph_left_32 = (def_j1(:, 2) .* b31) + (def_j1(:, 5) .* b32) + (def_j1(:, 8) .* b33);
                ph_left_33 = (def_j1(:, 3) .* b31) + (def_j1(:, 6) .* b32) + (def_j1(:, 9) .* b33);
 
                ph_left_def_j1 = [ph_left_11, ph_left_12, ph_left_13,...
                                  ph_left_21, ph_left_22, ph_left_23,...
                                  ph_left_31, ph_left_32, ph_left_33];
                clear def_j1





                ph_left_11 = (def_k1(:, 1) .* b11) + (def_k1(:, 4) .* b12) + (def_k1(:, 7) .* b13);
                ph_left_12 = (def_k1(:, 2) .* b11) + (def_k1(:, 5) .* b12) + (def_k1(:, 8) .* b13);
                ph_left_13 = (def_k1(:, 3) .* b11) + (def_k1(:, 6) .* b12) + (def_k1(:, 9) .* b13);
                ph_left_21 = (def_k1(:, 1) .* b21) + (def_k1(:, 4) .* b22) + (def_k1(:, 7) .* b23);
                ph_left_22 = (def_k1(:, 2) .* b21) + (def_k1(:, 5) .* b22) + (def_k1(:, 8) .* b23);
                ph_left_23 = (def_k1(:, 3) .* b21) + (def_k1(:, 6) .* b22) + (def_k1(:, 9) .* b23);
                ph_left_31 = (def_k1(:, 1) .* b31) + (def_k1(:, 4) .* b32) + (def_k1(:, 7) .* b33);
                ph_left_32 = (def_k1(:, 2) .* b31) + (def_k1(:, 5) .* b32) + (def_k1(:, 8) .* b33);
                ph_left_33 = (def_k1(:, 3) .* b31) + (def_k1(:, 6) .* b32) + (def_k1(:, 9) .* b33);
 
                ph_left_def_k1 = [ph_left_11, ph_left_12, ph_left_13,...
                                  ph_left_21, ph_left_22, ph_left_23,...
                                  ph_left_31, ph_left_32, ph_left_33];
                 clear def_k1

                ph_right_11 = (ph_left_def_j1(:, 1) .* b11) + (ph_left_def_j1(:, 2) .* b21) + (ph_left_def_j1(:, 3) .* b31);
                ph_right_12 = (ph_left_def_j1(:, 1) .* b12) + (ph_left_def_j1(:, 2) .* b22) + (ph_left_def_j1(:, 3) .* b32);
                ph_right_13 = (ph_left_def_j1(:, 1) .* b13) + (ph_left_def_j1(:, 2) .* b23) + (ph_left_def_j1(:, 3) .* b33);
                ph_right_21 = (ph_left_def_j1(:, 4) .* b11) + (ph_left_def_j1(:, 5) .* b21) + (ph_left_def_j1(:, 6) .* b31);
                ph_right_22 = (ph_left_def_j1(:, 4) .* b12) + (ph_left_def_j1(:, 5) .* b22) + (ph_left_def_j1(:, 6) .* b32);
                ph_right_23 = (ph_left_def_j1(:, 4) .* b13) + (ph_left_def_j1(:, 5) .* b23) + (ph_left_def_j1(:, 6) .* b33);
                ph_right_31 = (ph_left_def_j1(:, 7) .* b11) + (ph_left_def_j1(:, 8) .* b21) + (ph_left_def_j1(:, 9) .* b31);
                ph_right_32 = (ph_left_def_j1(:, 7) .* b12) + (ph_left_def_j1(:, 8) .* b22) + (ph_left_def_j1(:, 9) .* b32);
                ph_right_33 = (ph_left_def_j1(:, 7) .* b13) + (ph_left_def_j1(:, 8) .* b23) + (ph_left_def_j1(:, 9) .* b33);
 
                full_def_j_1 = [ph_right_11, ph_right_12, ph_right_13,...
                               ph_right_21, ph_right_22, ph_right_23,...
                               ph_right_31, ph_right_32, ph_right_33];
                
                 clear ph_left_def_j1




                ph_right_11 = (ph_left_def_k1(:, 1) .* b11) + (ph_left_def_k1(:, 2) .* b21) + (ph_left_def_k1(:, 3) .* b31);
                ph_right_12 = (ph_left_def_k1(:, 1) .* b12) + (ph_left_def_k1(:, 2) .* b22) + (ph_left_def_k1(:, 3) .* b32);
                ph_right_13 = (ph_left_def_k1(:, 1) .* b13) + (ph_left_def_k1(:, 2) .* b23) + (ph_left_def_k1(:, 3) .* b33);
                ph_right_21 = (ph_left_def_k1(:, 4) .* b11) + (ph_left_def_k1(:, 5) .* b21) + (ph_left_def_k1(:, 6) .* b31);
                ph_right_22 = (ph_left_def_k1(:, 4) .* b12) + (ph_left_def_k1(:, 5) .* b22) + (ph_left_def_k1(:, 6) .* b32);
                ph_right_23 = (ph_left_def_k1(:, 4) .* b13) + (ph_left_def_k1(:, 5) .* b23) + (ph_left_def_k1(:, 6) .* b33);
                ph_right_31 = (ph_left_def_k1(:, 7) .* b11) + (ph_left_def_k1(:, 8) .* b21) + (ph_left_def_k1(:, 9) .* b31);
                ph_right_32 = (ph_left_def_k1(:, 7) .* b12) + (ph_left_def_k1(:, 8) .* b22) + (ph_left_def_k1(:, 9) .* b32);
                ph_right_33 = (ph_left_def_k1(:, 7) .* b13) + (ph_left_def_k1(:, 8) .* b23) + (ph_left_def_k1(:, 9) .* b33);
 
                full_def_k_1 = [ph_right_11, ph_right_12, ph_right_13,...
                               ph_right_21, ph_right_22, ph_right_23,...
                               ph_right_31, ph_right_32, ph_right_33];
                
                clear ph_left_def_k1

                b11 = Ph_elem_edge(1, 1, m, :, 2);
                b11 = b11(:);
                b12 = Ph_elem_edge(1, 2, m, :, 2);
                b12 = b12(:);
                b13 = Ph_elem_edge(1, 3, m, :, 2);
                b13 = b13(:);
                b21 = Ph_elem_edge(2, 1, m, :, 2);
                b21 = b21(:);
                b22 = Ph_elem_edge(2, 2, m, :, 2);
                b22 = b22(:);
                b23 = Ph_elem_edge(2, 3, m, :, 2);
                b23 = b23(:);
                b31 = Ph_elem_edge(3, 1, m, :, 2);
                b31 = b31(:);
                b32 = Ph_elem_edge(3, 2, m, :, 2);
                b32 = b32(:);
                b33 = Ph_elem_edge(3, 3, m, :, 2);
                b33 = b33(:);

                ph_left_11 = (def_j2(:, 1) .* b11) + (def_j2(:, 4) .* b12) + (def_j2(:, 7) .* b13);
                ph_left_12 = (def_j2(:, 2) .* b11) + (def_j2(:, 5) .* b12) + (def_j2(:, 8) .* b13);
                ph_left_13 = (def_j2(:, 3) .* b11) + (def_j2(:, 6) .* b12) + (def_j2(:, 9) .* b13);
                ph_left_21 = (def_j2(:, 1) .* b21) + (def_j2(:, 4) .* b22) + (def_j2(:, 7) .* b23);
                ph_left_22 = (def_j2(:, 2) .* b21) + (def_j2(:, 5) .* b22) + (def_j2(:, 8) .* b23);
                ph_left_23 = (def_j2(:, 3) .* b21) + (def_j2(:, 6) .* b22) + (def_j2(:, 9) .* b23);
                ph_left_31 = (def_j2(:, 1) .* b31) + (def_j2(:, 4) .* b32) + (def_j2(:, 7) .* b33);
                ph_left_32 = (def_j2(:, 2) .* b31) + (def_j2(:, 5) .* b32) + (def_j2(:, 8) .* b33);
                ph_left_33 = (def_j2(:, 3) .* b31) + (def_j2(:, 6) .* b32) + (def_j2(:, 9) .* b33);
 
                ph_left_def_j2 = [ph_left_11, ph_left_12, ph_left_13,...
                                  ph_left_21, ph_left_22, ph_left_23,...
                                  ph_left_31, ph_left_32, ph_left_33];
                clear def_j1





                ph_left_11 = (def_k2(:, 1) .* b11) + (def_k2(:, 4) .* b12) + (def_k2(:, 7) .* b13);
                ph_left_12 = (def_k2(:, 2) .* b11) + (def_k2(:, 5) .* b12) + (def_k2(:, 8) .* b13);
                ph_left_13 = (def_k2(:, 3) .* b11) + (def_k2(:, 6) .* b12) + (def_k2(:, 9) .* b13);
                ph_left_21 = (def_k2(:, 1) .* b21) + (def_k2(:, 4) .* b22) + (def_k2(:, 7) .* b23);
                ph_left_22 = (def_k2(:, 2) .* b21) + (def_k2(:, 5) .* b22) + (def_k2(:, 8) .* b23);
                ph_left_23 = (def_k2(:, 3) .* b21) + (def_k2(:, 6) .* b22) + (def_k2(:, 9) .* b23);
                ph_left_31 = (def_k2(:, 1) .* b31) + (def_k2(:, 4) .* b32) + (def_k2(:, 7) .* b33);
                ph_left_32 = (def_k2(:, 2) .* b31) + (def_k2(:, 5) .* b32) + (def_k2(:, 8) .* b33);
                ph_left_33 = (def_k2(:, 3) .* b31) + (def_k2(:, 6) .* b32) + (def_k2(:, 9) .* b33);
 
                ph_left_def_k2 = [ph_left_11, ph_left_12, ph_left_13,...
                                  ph_left_21, ph_left_22, ph_left_23,...
                                  ph_left_31, ph_left_32, ph_left_33];
                 clear def_k1

                ph_right_11 = (ph_left_def_j2(:, 1) .* b11) + (ph_left_def_j2(:, 2) .* b21) + (ph_left_def_j2(:, 3) .* b31);
                ph_right_12 = (ph_left_def_j2(:, 1) .* b12) + (ph_left_def_j2(:, 2) .* b22) + (ph_left_def_j2(:, 3) .* b32);
                ph_right_13 = (ph_left_def_j2(:, 1) .* b13) + (ph_left_def_j2(:, 2) .* b23) + (ph_left_def_j2(:, 3) .* b33);
                ph_right_21 = (ph_left_def_j2(:, 4) .* b11) + (ph_left_def_j2(:, 5) .* b21) + (ph_left_def_j2(:, 6) .* b31);
                ph_right_22 = (ph_left_def_j2(:, 4) .* b12) + (ph_left_def_j2(:, 5) .* b22) + (ph_left_def_j2(:, 6) .* b32);
                ph_right_23 = (ph_left_def_j2(:, 4) .* b13) + (ph_left_def_j2(:, 5) .* b23) + (ph_left_def_j2(:, 6) .* b33);
                ph_right_31 = (ph_left_def_j2(:, 7) .* b11) + (ph_left_def_j2(:, 8) .* b21) + (ph_left_def_j2(:, 9) .* b31);
                ph_right_32 = (ph_left_def_j2(:, 7) .* b12) + (ph_left_def_j2(:, 8) .* b22) + (ph_left_def_j2(:, 9) .* b32);
                ph_right_33 = (ph_left_def_j2(:, 7) .* b13) + (ph_left_def_j2(:, 8) .* b23) + (ph_left_def_j2(:, 9) .* b33);
 
                full_def_j_2 = [ph_right_11, ph_right_12, ph_right_13,...
                               ph_right_21, ph_right_22, ph_right_23,...
                               ph_right_31, ph_right_32, ph_right_33];
                
                 clear ph_left_def_j1




                ph_right_11 = (ph_left_def_k2(:, 1) .* b11) + (ph_left_def_k2(:, 2) .* b21) + (ph_left_def_k2(:, 3) .* b31);
                ph_right_12 = (ph_left_def_k2(:, 1) .* b12) + (ph_left_def_k2(:, 2) .* b22) + (ph_left_def_k2(:, 3) .* b32);
                ph_right_13 = (ph_left_def_k2(:, 1) .* b13) + (ph_left_def_k2(:, 2) .* b23) + (ph_left_def_k2(:, 3) .* b33);
                ph_right_21 = (ph_left_def_k2(:, 4) .* b11) + (ph_left_def_k2(:, 5) .* b21) + (ph_left_def_k2(:, 6) .* b31);
                ph_right_22 = (ph_left_def_k2(:, 4) .* b12) + (ph_left_def_k2(:, 5) .* b22) + (ph_left_def_k2(:, 6) .* b32);
                ph_right_23 = (ph_left_def_k2(:, 4) .* b13) + (ph_left_def_k2(:, 5) .* b23) + (ph_left_def_k2(:, 6) .* b33);
                ph_right_31 = (ph_left_def_k2(:, 7) .* b11) + (ph_left_def_k2(:, 8) .* b21) + (ph_left_def_k2(:, 9) .* b31);
                ph_right_32 = (ph_left_def_k2(:, 7) .* b12) + (ph_left_def_k2(:, 8) .* b22) + (ph_left_def_k2(:, 9) .* b32);
                ph_right_33 = (ph_left_def_k2(:, 7) .* b13) + (ph_left_def_k2(:, 8) .* b23) + (ph_left_def_k2(:, 9) .* b33);
 
                full_def_k_2 = [ph_right_11, ph_right_12, ph_right_13,...
                               ph_right_21, ph_right_22, ph_right_23,...
                               ph_right_31, ph_right_32, ph_right_33];
                
                clear ph_left_def_k1

                % full_def_j_1 full_def_j_2 full_def_k_1 full_def_k_2
                t1normal1 = t_conormal_edge(1, 1, m, :, 1);
                t1normal1 = t1normal1(:);
                t1normal2 = t_conormal_edge(1, 2, m, :, 1);
                t1normal2 = t1normal2(:);
                t1normal3 = t_conormal_edge(1, 3, m, :, 1);
                t1normal3 = t1normal3(:);
                t2normal1 = t_conormal_edge(1, 1, m, :, 2);
                t2normal1 = t2normal1(:);
                t2normal2 = t_conormal_edge(1, 2, m, :, 2);
                t2normal2 = t2normal2(:);
                t2normal3 = t_conormal_edge(1, 3, m, :, 2);
                t2normal3 = t2normal3(:);
                edgetangent11 = edge_tangent(1, 1, m, :);
                edgetangent11 = edgetangent11(:);
                edgetangent12 = edge_tangent(1, 2, m, :);
                edgetangent12 = edgetangent12(:);
                edgetangent13 = edge_tangent(1, 3, m, :);
                edgetangent13 = edgetangent13(:);
                hedge2 = h_edge(:, m, :);
                hedge2 = hedge2(:);

                dgamj1 = edgetangent11 .* sum(full_def_j_1(:, 1:3).* [t1normal1, t1normal2, t1normal3], 2) +...
                         edgetangent12 .* sum(full_def_j_1(:, 4:6).* [t1normal1, t1normal2, t1normal3], 2) +...
                         edgetangent13 .* sum(full_def_j_1(:, 7:9).* [t1normal1, t1normal2, t1normal3], 2);

                dgamj2 = edgetangent11 .* sum(full_def_j_2(:, 1:3).* [t2normal1, t2normal2, t2normal3], 2) +...
                         edgetangent12 .* sum(full_def_j_2(:, 4:6).* [t2normal1, t2normal2, t2normal3], 2) +...
                         edgetangent13 .* sum(full_def_j_2(:, 7:9).* [t2normal1, t2normal2, t2normal3], 2);
                dgamk1 = edgetangent11 .* sum(full_def_k_1(:, 1:3).* [t1normal1, t1normal2, t1normal3], 2) +...
                         edgetangent12 .* sum(full_def_k_1(:, 4:6).* [t1normal1, t1normal2, t1normal3], 2) +...
                         edgetangent13 .* sum(full_def_k_1(:, 7:9).* [t1normal1, t1normal2, t1normal3], 2);
                dgamk2 = edgetangent11 .* sum(full_def_k_2(:, 1:3).* [t2normal1, t2normal2, t2normal3], 2) +...
                         edgetangent12 .* sum(full_def_k_2(:, 4:6).* [t2normal1, t2normal2, t2normal3], 2) +...
                         edgetangent13 .* sum(full_def_k_2(:, 7:9).* [t2normal1, t2normal2, t2normal3], 2);

                A1jk  = A1jk  + edgequad.qwts(m) .* jds1 .* kds1 .*...
                              (-(dgamj1 .* u1k + dgamk1 .* u1j) .*...
                              hedge2 ./ 2 + beta .* u1j .* u1k);
                

                A2jk = A2jk + edgequad.qwts(m) .* jds2 .* kds2 .*...
                              (-(dgamj2 .* u2k + dgamk2 .* u2j) .*...
                              hedge2 ./ 2 + beta .* u2j .* u2k);


                A3jk = A3jk + edgequad.qwts(m) .* jds1 .* kds2 .*...
                              ((dgamj1 .* u2k + dgamk2 .* u1j) .*...
                              hedge2 ./ 2 - beta .*  u1j .* u2k);
        

                A4jk = A4jk + edgequad.qwts(m) .* kds1 .* jds2 .*...
                              ((dgamj2 .* u1k + dgamk1 .* u2j) .*...
                              hedge2 ./ 2 - beta .* u2j .* u1k);
            end
            AAjindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 1), j);
            AAkindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 1), k);
            AAvec(iter*numedge+1:(iter+1)*numedge) = A1jk;
            iter = iter + 1;
            AAjindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 2), j);
            AAkindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 2), k);
            AAvec(iter*numedge+1:(iter+1)*numedge) = A2jk;
            iter = iter + 1;
            AAjindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 1), j);
            AAkindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 2), k);
            AAvec(iter*numedge+1:(iter+1)*numedge) = A3jk;
            iter = iter + 1; 
            AAjindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 2), j);
            AAkindex(iter*numedge+1:(iter+1)*numedge) = elemdof(edge2elem(:, 1), k);
            AAvec(iter*numedge+1:(iter+1)*numedge) = A4jk;
            iter = iter + 1;

        end        
    end


toc
A = sparse([Ajindex; AAjindex], [Akindex; AAkindex], ...
         [Avec; AAvec], 3*numedge+3*numel, 3*numedge+3*numel);


% %%


A = (A + A') / 2;
M = (M + M') / 2;
%  
% left_top_part = 2 * A +  M;
% first_part_schur_rhs = left_top_part\rhs;
% %So we have the schur right hand side
% schur_rhs = B' * first_part_schur_rhs - grhs;
% schur_matrix = B' /(2 * A + 1 * M) * B;
% pvec = pcg(schur_matrix, schur_rhs, 10^(-7), 5000);
% %uvec = pcg(2 * A + h * M, rhs - B * pvec, 10^(-7), 5000);
% uvec = (2 * A + 1 * M)\(rhs - B * pvec);

upvec = [2*A+M B; B' sparse(3*numel,3*numel)]\[rhs; grhs];
upvec = full(upvec);
uvec = upvec(1:3*numedge+3*numel);
[L, U] = ilu(B'*B);
pvec = pcg(B'*B, B'*rhs-B'*(2*A+M)*uvec, 10^(-7), 5000, L, U);


vel_l2err = 0;
grad_error_vel = 0;
l2_press_err = 0;
for m=1:recquad.nqpts
    
    
    for i=1:numel
        p = noden(i, :, m);
        proj_point_surface = surfacedata.project(p);
        press_val = -psurf(proj_point_surface);
        uval = usurf_perp_el(proj_point_surface, aa)';
        nuh = mycross(grad_ak(:, 1, m, i)', grad_ak(:, 2, m, i)');
        nuh = nuh / norm(nuh, 2);
        nu = surfacedata.unitoutnormal(proj_point_surface);
        nudotnuh = nu * nuh';
        piola_transform = eye(3) - nu' * nuh / nudotnuh;
        puval = piola_transform * uval;
        uhloc = zeros(3, 1);
        grad_val = jacob_usurf_perp_el(proj_point_surface, aa);
        grad_val = reshape(grad_val, 3, 3);
        grad_val = piola_transform * grad_val';
        grad_loc = zeros(3, 3);
        Ph = eye(3) - nuh' * nuh;
        phloc = 0;
        for j=1:12
            jdofind = elemdof(i, j);
            jdofsign = elemdofsign(i, j);
            uhloc = uhloc + uvec(jdofind) * jdofsign * grad_ak(:, :, m, i) *...
                    bdm_2basis_orsan(recquad.qpts(m, :), j)' / jac_det_elem(:, m, i); 
            jacob1 = hessian_matrix(:, :, m, i, 1);
            jacob2 = hessian_matrix(:, :, m, i, 2);
            jacob3 = hessian_matrix(:, :, m, i, 3);
            a111 = jacob1(1, 1);
            a112 = jacob1(1, 2);
            a121 = a112;
            a122 = jacob1(2, 2);
            a211 = jacob2(1, 1);
            a212 = jacob2(1, 2);
            a221 = a212;
            a222 = jacob2(2, 2);
            a311 = jacob3(1, 1);
            a312 = jacob3(1, 2);
            a321 = a312;
            a322 = jacob3(2, 2);
            a11 = grad_ak(1, 1, m, i);
            a12 = grad_ak(1, 2, m, i);
            a21 = grad_ak(2, 1, m, i);
            a22 = grad_ak(2, 2, m, i);
            a31 = grad_ak(3, 1, m, i);
            a32 = grad_ak(3, 2, m, i);
% 
            jacob_matrixj = [bdm_2basis_orsan(recquad.qpts(m, :), j)*jacob1;
                             bdm_2basis_orsan(recquad.qpts(m, :), j)*jacob2;
                             bdm_2basis_orsan(recquad.qpts(m, :), j)*jacob3]/jac_det_elem(:, m, i); 

            pseudo_inv_comp = (grad_ak(:, :, m, i)'*grad_ak(:, :, m, i));             
            gradient_partj = (grad_ak(:, :, m, i) * bdm2grad_basis_orsan(recquad.qpts(m, :), j))/jac_det_elem(:, m, i); 
                            det_grad_part1 = 2*a111*a11*a22^2 + 2*a121*a12*a21^2 + 2*a111*a11*a32^2 + 2*a211*a12^2*a21 + 2*a121*a12*a31^2 + 2*a221*a11^2*a22 + 2*a211*a21*a32^2 + 2*a311*a12^2*a31 + 2*a221*a22*a31^2 + 2*a321*a11^2*a32 + 2*a311*a22^2*a31 + 2*a321*a21^2*a32 - 2*a111*a12*a21*a22 - 2*a121*a11*a21*a22 - 2*a211*a11*a12*a22 - 2*a221*a11*a12*a21 - 2*a111*a12*a31*a32 - 2*a121*a11*a31*a32 - 2*a311*a11*a12*a32 - 2*a321*a11*a12*a31 - 2*a211*a22*a31*a32 - 2*a221*a21*a31*a32 - 2*a311*a21*a22*a32 - 2*a321*a21*a22*a31;
            det_grad_part2 = 2*a112*a11*a22^2 + 2*a122*a12*a21^2 + 2*a112*a11*a32^2 + 2*a212*a12^2*a21 + 2*a122*a12*a31^2 + 2*a222*a11^2*a22 + 2*a212*a21*a32^2 + 2*a312*a12^2*a31 + 2*a222*a22*a31^2 + 2*a322*a11^2*a32 + 2*a312*a22^2*a31 + 2*a322*a21^2*a32 - 2*a112*a12*a21*a22 - 2*a122*a11*a21*a22 - 2*a212*a11*a12*a22 - 2*a222*a11*a12*a21 - 2*a112*a12*a31*a32 - 2*a122*a11*a31*a32 - 2*a312*a11*a12*a32 - 2*a322*a11*a12*a31 - 2*a212*a22*a31*a32 - 2*a222*a21*a31*a32 - 2*a312*a21*a22*a32 - 2*a322*a21*a22*a31;

            det_grad1 = [det_grad_part1 det_grad_part2];
            det_part1 = -1/2*grad_ak(:, :, m, i) * bdm_2basis_orsan(recquad.qpts(m, :), j)'*det_grad1/jac_det_elem(:, m, i)^3;
            grad_j = (jacob_matrixj+det_part1+gradient_partj) / pseudo_inv_comp * grad_ak(:, :, m, i)';
            grad_loc = grad_loc + uvec(jdofind) * jdofsign * grad_j;
        end
        for t=1:3
            tdofind = press_dof(i, t);
            phloc = phloc + pvec(tdofind) *...
                    p1basis_orsan(recquad.qpts(m, :), t)/ jac_det_elem(:, m, i);
        end
        uloc = puval;
        ploc = press_val;
        grad_val_loc = grad_val;
        l2_press_err = l2_press_err + ...
                       (ploc - phloc)^2 * jac_det_elem(:, m, i) * recquad.qwts(m);
        vel_l2err = vel_l2err + ...
                (uloc - uhloc)' * (uloc - uhloc) * jac_det_elem(:, m, i) * recquad.qwts(m);
        grad_error_vel = grad_error_vel + ...
                 sum(dot(Ph * (grad_loc - grad_val_loc) * Ph, Ph * (grad_loc - grad_val_loc) * Ph), 2) * jac_det_elem(:, m, i) * recquad.qwts(m);
        
    end
end

%
finall2err_epsh = sqrt(vel_l2err)
final_grad_error_vel = sqrt(grad_error_vel)
final_l2err_press = sqrt(l2_press_err)