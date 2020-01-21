function[Temp_2,Temp_3,Temp_4,Temp_87,Temp_88,Temp_89] = FDM_TempModel_GSA(matrix,x,ProcessParam,MaterialParam)

%Process Variables
T_L = ProcessParam(1); %Extrusion temperature (ºC)
T_E = ProcessParam(2); %Temperature of the envelope (ºC)
v = ProcessParam(3); % printer speed/Velocity of the extrusion head (m/sec)
%Filament dimensions
Wid = ProcessParam(4); % Layer Thickness (meters)
Ht = ProcessParam(5); % Layer Height (meters)
L = ProcessParam(6); % Length of the filament (meters)

ro = MaterialParam(1); % Density (kg/m^3)
C = MaterialParam(2); % Specific heat (J/kg.K)
h_conv = MaterialParam(3); % Heat transfer coefficient (loss of heat by natural convection)
% h_cond = MaterialParam(:,4); % filament and adjacent filament & support


lam = MaterialParam(4); % filament and adjacent filament & support


% N_Samples = size(ro,1);

                
%____________________________________ STEP 1 ____________________________________
%Definition of the vector that contains the number of total filaments in each layer
matrix_lin = size(matrix,1);
matrix_col = size(matrix,2);
vector = zeros(matrix_lin,2);
ctr = 0; % # of layers
for i = matrix_lin:-1:1
    ctr = ctr + 1;
    for j = 1:matrix_col
        if matrix(i,j) ~= 0
            vector(ctr,1) = vector(ctr,1) + 1;
        end
    end
end

%Number of layers
m = length(vector(:,1));
%Number of filaments
n = 0;
for j = 1:m
    if m == 1
        n = vector(1,1);
    else
        if vector(j,2) <= 1
            n = n + vector(j,1);
        end
    end
end
%____________________________________ STEP 2 ____________________________________
%Computation variables
dt = .1; %Step time
temp_extra = 10; %Additional time computation after construction of the part
erro = 1e-3; %Convergence error
% FinalTimeFirstFila = .1; % delay time between filamentss
%____________________________________ STEP 3 ____________________________________
%Definition of the size of the variables
% ii = 5; # of contacts with adjacent filaments and support
h = zeros(1,5); lambda = zeros(1,5); a = zeros(n,5); T = zeros (n,5);
vec_b = zeros(n,5); vec_Q = zeros(n,5); b = zeros(1,n); Q = zeros(1,n);
T_begin = zeros(1,n); dif = zeros(1,n); Biot = zeros(1,n); save_T = zeros(1,n);
old_T = zeros(1,n); save_lim = zeros(1,n); viz = zeros(11,n);
%____________________________________ STEP 4 ____________________________________
%Process Variables
for lin = 1:n %Temperature of support (ºC)
    T(lin,5) = T_E;
end

% % For Cylinder
% area = pi * (Wid/2)^2; %Area of the cross section of filament (meters^2) -
% per = pi * Wid; %Perimeter of the cross section of filament (meters)

% For Ellipse
aa=Wid/2; bb=Ht/2;
hh = (aa-bb)^2/(aa+bb)^2;
area = pi * (aa*bb); %Ellipse:Area of the cross section of filament (meters^2)
per = pi * (aa+bb)*(1+3*hh/(10+sqrt(4-3*hh))); %Perimeter of the cross section of filament (meters)

vol = area*L; %Volume of the filament
A_p = per*L; %Superficial area of the filament

% Material Properties
%Thermal conductivity (W/m.K)
% conductivity(1) = 0.1768; % material A
% conductivity(1) = 0.15; % material A
% conductivity(2) = 0.5; % material B
conductivity = 0.15; % material A
% %Density (kg/m^3)
% ro(1) = 1050; % material A
% ro(1) = 1040; % material A
% ro(2) = 1500; % material B
% %Specific heat (J/kg.K)
% C(1) = 2019.7; % material A
% C(1) = 1290; % material A
% C(2) = 2500.7; % material B
%____________________________________ STEP 5 ____________________________________
% Heat transfer coefficient (loss of heat by natural convection)
% h_conv = 86;
%Thermal contact conductances between
h(1,1) = 86; % filament and left adjacent filament
h(1,2) = 86; % filament and down adjacent filament
h(1,3) = 86; % filament and right adjacent filament
h(1,4) = 86; % filament and top adjacent filament
h(1,5) = 86; % filament and support
% h(1,5) = h_cond;
%Fraction of perimeter contact between
lambda(1,1) = 0.1; % filament and left adjacent filament
lambda(1,2) = 0.1; % filament and down adjacent filament
lambda(1,3) = 0.2; % filament and right adjacent filament
lambda(1,4) = 0.1; % filament and top adjacent filament
lambda(1,5) = 0.1; % filament and support
lambda(1,1:5)=lam;
%____________________________________ STEP 6 ____________________________________
%Definition of the parameters influenced by the contacts
for col = 1:5
    for lin = 1:n
        vec_b(lin,col) = h(1,col)*lambda(1,col);
        vec_Q(lin,col) = vec_b(lin,col)*T(lin,col);
    end
end
% %____________________________________ STEP 7 ____________________________________
% %Definition of the parameters influenced by the material properties
% ctr = 0;
% number_filament = 0;
% for i = matrix_lin:-1:1
%     ctr = ctr + 1;
%     % mod(N,2)==0 % iseven
%     % mod(N,2)==1 % isodd
%     %     if isodd(contar) == 1
%     if mod(ctr,2) == 1
%         for j = 1:matrix_col
%             if matrix(i,j) ~= 0
%                 number_filament = number_filament + 1;
%                 scalar(number_filament) = -per/(ro(matrix(i,j))*area*C(matrix(i,j)));
%                 esc(number_filament) = h_conv/(ro(matrix(i,j))*L*C(matrix(i,j)));
%                 kt(number_filament) = conductivity(matrix(i,j));
%             end
%         end
%     else
%         for j = matrix_col:-1:1
%             if matrix(i,j) ~= 0
%                 number_filament = number_filament + 1;
%                 scalar(number_filament) = -per/(ro(matrix(i,j))*area*C(matrix(i,j)));
%                 esc(number_filament) = h_conv/(ro(matrix(i,j))*L*C(matrix(i,j)));
%                 kt(number_filament) = conductivity(matrix(i,j));
%             end
%         end
%     end
% end


scalar = -per/(ro*area*C);
kt = conductivity;

limite=zeros(n+2,2); % BK

for kk = 1:1:(n+2)
    if mod(kk,2) == 1
        limite(kk,1) = (kk*L-x)/v;
        limite(kk,2) = (kk*L+x)/v;
    else
        limite(kk,1) = limite(kk-1,2);
        limite(kk,2) = ((kk+1)*L-x)/v;
    end
end

% temp=zeros(); % BK
for road = 1:n
    line = 0;
    for i = 0:dt:limite(n,2)
        line = line + 1;
        temp(line,road) = T_L;
    end
end

% abcissa=zeros((limite(end,2)+temp_extra)/dt,1); % BK
for layer = 1:m
    if layer == 1
        for num = 1:vector(layer,1)
            if num == 1
                %____________________________________ STEP 9 ____________________________________
                a(num,5) = 1; %Activation of the contact with support
                %____________________________________ STEP 10 ____________________________________
                %Definition of the variables b and Q defined in equation Eq. 7
                b(num) = h_conv*(1-lambda*a(num,:)') + vec_b(num,:)*a(num,:)';
                Q(num) = (h_conv*(1-lambda*a(num,:)')*T_E + ...
                    vec_Q(num,:)*a(num,:)')/b(num);
                
                %____________________________________ STEP 11 ____________________________________
                p = 0;
                for t = 0:dt:limite(num,1)
                    p = p+1; abcissa(p) = t;
                end
                %____________________________________ STEP 12 ____________________________________
                %Computation of the temperatures of the first filament
                %                 for t = (limite(num,1)+dt):dt:limite_final
                %                 for t = (limite(num,1)+dt):dt:FinalTimeFirstFila % BK
                for t = (limite(num,1)+dt):dt:limite(1,2) % BK
                    p = p+1; abcissa(p) = t;
                    temp(p,num)=(T_L-Q(num))*exp(scalar(num)*b(num)* ...
                        (t-limite(num,1))) + Q(num);
                end
                %Saving the last temperature of the period time of cooling down
                T_begin(num) = temp(p,num);
                %____________________________________ STEP 13 ____________________________________
                %Verification of the value of Biot Number
                Biot(num) = (vol/A_p)*(b(num)/kt(num));
                if Biot(num)>=0.1
                    'WARNING! We cannot use a Lumped System';
                end
                
                %____________________________________ STEP 14 ____________________________________
            else
                %Activation of the contacts
                a(num-1,3) = 1; a(num,1) = 1; a(num,5) = 1;
                %____________________________________ STEP 15 ____________________________________
                %Up-dating of the variable b
                for j = 1:num
                    b(j) = h_conv*(1-lambda*a(j,:)') + vec_b(j,:)*a(j,:)';
                end
                
                %____________________________________ STEP 16 ____________________________________
                if m == 1
                    if num == vector(layer,1)
                        limite_final = limite(num,2) + temp_extra;
                    else
                        limite_final = limite(num,2);
                    end
                else
                    limite_final = limite(num,2);
                end
                
                %                 for t = (FinalTimeFirstFila+dt):dt:limite_final % BK
                for t = (limite(num,1)+dt):dt:limite_final
                    p = p+1; abcissa(p) = t;
                    last = p-1;
                    for j = 1:num
                        save_T(j) = temp(last,j);
                    end
                    %____________________________________ STEP 17 ____________________________________
                    %Iterative process
                    for q = 1:100000
                        %Saving contacts and temperatures of adjacent filaments
                        for j = 1:num
                            if j == 1
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                            end
                            if j > 1 && j < num
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                            end
                            if j == num
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                            end
                            for k = 1:5
                                if T(j,k) ~= 0 && k ~= 5
                                    vec_Q(j,k) = vec_b(j,k)*T(j,k);
                                end
                            end
                            %Up-dating of the variable Q
                            Q(j) = (h_conv*(1-lambda*a(j,:)')*T_E + ...
                                vec_Q(j,:)*a(j,:)')/b(j);
                            old_T(j) = save_T(j);
                        end
                        %Computation of the temperatures
                        if num == 2
                            save_T(num-1) = (T_begin(num-1)-Q(num-1))* ...
                                exp(scalar*b(num-1)*(t-limite(num,1)))+Q(num-1);
                            save_T(num) = (T_L-Q(num))* ...
                                exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            save_lim(1,1) = limite(num,1);
                            save_lim(1,2) = limite(num,1);
                        else
                            for j=1:num-2
                                save_T(j) = (T_begin(j)-Q(j))*exp(scalar*b(j)* ...
                                    (t-limite(j,1)))+Q(j);
                            end
                            save_T(num-1) = (T_begin(num-1)-Q(num-1))* ...
                                exp(scalar*b(num-1)*(t-limite(num,1)))+Q(num-1);
                            save_T(num) = (T_L-Q(num))* ...
                                exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            save_lim(1,num-1) = limite(num,1);
                            save_lim(1,num) = limite(num,1);
                        end
                        for j = 1:num
                            dif(j) = abs(save_T(j)-old_T(j));
                        end
                        try_ = 1;
                        stop = 0;
                        for j = 1:num
                            if dif(try_) < erro
                                try_ = try_+1;
                            end
                            if try_ == num+1
                                stop = 1;
                            end
                        end
                        if stop == 1
                            for j = 1:num
                                temp(p,j) = save_T(j);
                            end
                            break;
                        end
                    end
                end
                T_begin(num) = temp(p,num);
                %End of iterative process
                %____________________________________ STEP 18 ____________________________________
                %Verification of the Biot Number
                for j=1:num
                    Biot(j) = (vol/A_p)*(b(j)/kt);
                    if Biot(j)>=0.1
                        'WARNING! We can not use a Lumped System';
                        j;
                        Biot(j);
                    end
                end
            end
        end
        % =================================================================== %
        % For the remaining layers -> BK
    else
        jj=1;
        for num = sum(vector(1:layer-1,1))+1:sum(vector(1:layer-1,1))+vector(layer,1)
            
            % =================================================================== %
            if mod(layer,2)==0 % for even layers
                
                % =================================================================== %
                % if the filament is the first filament on that layer
                if num == sum(vector(1:layer-1,1))+1
                    %____________________________________ STEP 9 ____________________________________
                    a(num,2) = 1; %Activation of the contact down filament
                    a(num-1,4) = 1; %Activation of the contact top filament
                    
                    % rest of the filaments on that layer
                else % if the filament is not the first filament on that layer
                    jj=jj+1;
                    
                    %____________________________________ STEP 14 ____________________________________
                    %Activation of the contacts (horizontal)
                    a(num-1,1) = 1; a(num,3) = 1;
                    %Activation of the contacts (vertical)
                    a(num,2) = 1; a(num-(2*jj-1),4) = 1;
                end
                
                %____________________________________ STEP 15 ____________________________________
                %Up-dating of the variable b
                for j = 1:num
                    b(j) = h_conv*(1-lambda*a(j,:)') + vec_b(j,:)*a(j,:)';
                end
                
                %____________________________________ STEP 16 ____________________________________
                if m == layer
                    if num == sum(vector(:,1))
                        limite_final = limite(num,2) + temp_extra;
                    else
                        limite_final = limite(num,2);
                    end
                else
                    limite_final = limite(num,2);
                end
                
                % check is array element is empty or not
                chk = size(save_T,2);
                
                for t = (limite(num,1)+dt):dt:limite_final
                    p = p+1; abcissa(p) = t;
                    last = p-1;
                    for j = 1:num
                        save_T(j) = temp(last,j);
                    end
                    
                    %____________________________________ STEP 17 ____________________________________
                    %Iterative process
                    for q = 1:10000
                        %Saving contacts and temperatures of adjacent filaments
                        for j = 1:num
                            if j == 1
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                if chk >= j+2*(sum(vector(1,1))-j)+1
                                    T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                    viz(4,j) = j+2*(sum(vector(1,1))-j);
                                end
                            end
                            if j > 1 && j < sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                if chk >= j+2*(sum(vector(1,1))-j)+1
                                    T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                    viz(4,j) = j+2*(sum(vector(1,1))-j);
                                end
                            end
                            if j == sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                if chk >= j+1
                                    T(j,4) = save_T(j+1);
                                    viz(4,j) = j+1;
                                end
                            end
                            
                            % loop over all even layers except first one
                            for kk=2:2:m
                                % even layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                end
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                            end
                            
                            % loop over all odd layers except first one
                            for kk=3:2:m
                                % odd layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                    T(j,4) = save_T(j+(2*(j-(sum(vector(1:kk-1,1))+1))));
                                    viz(4,j) = j+(2*(j-(sum(vector(1:kk-1,1))+1)));
                                end
                                % odd layers' inner filaments
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    T(j,4) = save_T(j+(2*(j-(sum(vector(1:kk,1))+2))));
                                    viz(4,j) = j+(2*(j-(sum(vector(1:kk-1,1))+2)));
                                end
                                % odd layers' last filament
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    T(j,4) = save_T(j+1);
                                    viz(4,j) = j+1;
                                end
                            end
                            
                            for k = 1:5
                                if T(j,k) ~= 0 && k ~= 5
                                    vec_Q(j,k) = vec_b(j,k)*T(j,k);
                                end
                            end
                            %Up-dating of the variable Q
                            Q(j) = (h_conv*(1-lambda*a(j,:)')*T_E + ...
                                vec_Q(j,:)*a(j,:)')/b(j);
                            old_T(j) = save_T(j);
                        end
                        
                        %Computation of the temperatures
                        for j=1:num-2
                            save_T(j) = (T_begin(j)-Q(j))*exp(scalar*b(j)* ...
                                (t-limite(j,1)))+Q(j);
                        end
                        save_T(num-1) = (T_begin(num-1)-Q(num-1))*exp(scalar*b(num-1)* ...
                            (t-limite(num,1)))+Q(num-1);
                        save_T(num) = (T_L-Q(num))*exp(scalar*b(num)*(t-limite(num,1)))+Q(num);
                        save_lim(1,num-1) = limite(num,1);
                        save_lim(1,num) = limite(num,1);
                        
                        for j = 1:num
                            dif(j) = abs(save_T(j)-old_T(j));
                        end
                        try_ = 1;
                        stop = 0;
                        for j = 1:num
                            if dif(try_) < erro
                                try_ = try_+1;
                            end
                            if try_ == num+1
                                stop = 1;
                            end
                        end
                        if stop == 1
                            for j = 1:num
                                temp(p,j) = save_T(j);
                            end
                            break;
                        end
                    end
                end
                T_begin(num) = temp(p,num);
                %End of iterative process
                %____________________________________ STEP 18 ____________________________________
                %Verification of the Biot Number
                for j=1:num
                    Biot(j) = (vol/A_p)*(b(j)/kt);
                    if Biot(j)>=0.1
                        'WARNING! We can not use a Lumped System';
                        j;
                        Biot(j);
                    end
                end
                
                % =================================================================== %
            else % odd layers
                % =================================================================== %
                % if the filament is the first filament on that layer
                if num == sum(vector(1:layer-1,1))+1
                    %____________________________________ STEP 9 ____________________________________
                    a(num,2) = 1; %Activation of the contact down filament
                    a(num-1,4) = 1; %Activation of the contact top filament
                    
                    % rest of the filaments on that layer
                else % if the filament is not the first filament on that layer
                    jj=jj+1;
                    
                    %____________________________________ STEP 14 ____________________________________
                    %Activation of the contacts (horizontal)
                    a(num-1,3) = 1; a(num,1) = 1;
                    %Activation of the contacts (vertical)
                    a(num,2) = 1; a(num-(2*jj-1),4) = 1;
                end
                
                %____________________________________ STEP 15 ____________________________________
                %Up-dating of the variable b
                for j = 1:num
                    b(j) = h_conv*(1-lambda*a(j,:)') + vec_b(j,:)*a(j,:)';
                end
                
                %____________________________________ STEP 16 ____________________________________
                if m == layer
                    if num == sum(vector(:,1))
                        limite_final = limite(num,2) + temp_extra;
                    else
                        limite_final = limite(num,2);
                    end
                else
                    limite_final = limite(num,2);
                end
                
                % check is array element is empty or not
                chk = size(save_T,2);
                
                for t = (limite(num,1)+dt):dt:limite_final
                    p = p+1; abcissa(p) = t;
                    last = p-1;
                    for j = 1:num
                        save_T(j) = temp(last,j);
                    end
                    
                    %____________________________________ STEP 17 ____________________________________
                    %Iterative process
                    for q = 1:100000
                        %Saving contacts and temperatures of adjacent filaments
                        for j = 1:num
                            %Saving contacts and temperatures of adjacent filaments
                            if j == 1
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                viz(4,j) = j+2*(sum(vector(1,1))-j);
                            end
                            if j > 1 && j < sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                viz(4,j) = j+2*(sum(vector(1,1))-j);
                            end
                            if j == sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,4) = save_T(j+1);
                                viz(4,j) = j+1;
                            end
                            
                            % loop over all even layers except first one
                            for kk=2:2:m
                                % even layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                    T(j,4) = save_T(j+(2*((sum(vector(1:kk-1,1))-j))));
                                    viz(4,j) = j+(2*((sum(vector(1:kk-1,1))-j)));
                                end
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    %                                     j+(2*((sum(vector(1:kk,1))-j))+1)
                                    T(j,4) = save_T(j+(2*((sum(vector(1:kk-1,1))-j)))+1);
                                    viz(4,j) = j+(2*((sum(vector(1:kk-1,1))-j)));
                                end
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    T(j,4) = save_T(j+1);
                                    viz(4,j) = j+1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                end
                            end
                            
                            % loop over all odd layers except first one
                            for kk=3:2:m
                                % odd layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                end
                                % odd layers' inner filaments
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                                % odd layers' last filament
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                            end
                            
                            for k = 1:5
                                if T(j,k) ~= 0 && k ~= 5
                                    vec_Q(j,k) = vec_b(j,k)*T(j,k);
                                end
                            end
                            %Up-dating of the variable Q
                            Q(j) = (h_conv*(1-lambda*a(j,:)')*T_E + ...
                                vec_Q(j,:)*a(j,:)')/b(j);
                            old_T(j) = save_T(j);
                        end
                        %Computation of the temperatures
                        for j=1:num-2
                            save_T(j) = (T_begin(j)-Q(j))*exp(scalar*b(j)* ...
                                (t-limite(j,1)))+Q(j);
                        end
                        save_T(num-1) = (T_begin(num-1)-Q(num-1))*exp(scalar*b(num-1)* ...
                            (t-limite(num,1)))+Q(num-1);
                        save_T(num) = (T_L-Q(num))*exp(scalar*b(num)*(t-limite(num,1)))+Q(num);
                        save_lim(1,num-1) = limite(num,1);
                        save_lim(1,num) = limite(num,1);
                        
                        for j = 1:num
                            dif(j) = abs(save_T(j)-old_T(j));
                        end
                        try_ = 1;
                        stop = 0;
                        for j = 1:num
                            if dif(try_) < erro
                                try_ = try_+1;
                            end
                            if try_ == num+1
                                stop = 1;
                            end
                        end
                        if stop == 1
                            for j = 1:num
                                temp(p,j) = save_T(j);
                            end
                            break;
                        end
                    end
                end
                
                T_begin(num) = temp(p,num);
                %End of iterative process
                %____________________________________ STEP 18 ____________________________________
                %Verification of the Biot Number
                for j=1:num
                    Biot(j) = (vol/A_p)*(b(j)/kt);
                    if Biot(j)>=0.1
                        'WARNING! We can not use a Lumped System';
                        j;
                        Biot(j);
                    end
                end
                
            end
        end
    end
end

% only output selected Temp.
Temp_2 = temp(:,2)';
Temp_3 = temp(:,3)';
Temp_4 = temp(:,4)';


Temp_87 = temp(:,87)';
Temp_88 = temp(:,88)';
Temp_89 = temp(:,89)';


% [~,idx]=min(abs(abcissa-350));
% time=abcissa(1:idx);
% temp=temp(1:idx,:);
% % figure3: Plot for Temperatures of each interface
% figure
% set(groot,'DefaultAxesColorOrder',jet(n),...
%     'DefaultAxesLineStyleOrder','-|--|:|-.');
% hold on;
% for ii = 1:n
%     plot(time,temp(:,ii),'linewidth',1);
% end
% xlabel('Time (s)', 'FontName', 'Times New Roman', ...
%     'FontSize',16,'Color','k', 'Interpreter', 'tex')
% ylabel('Temperature ({\circ}C)', 'FontName', 'Times New Roman', ...
%     'FontSize',16,'Color','k', 'Interpreter', 'tex')
% set(gcf, 'PaperPosition', [0 0 4 6]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [4 6]); %Set the paper to have width 5 and height 5.
% % h = gca;  % Handle to currently active axes
% set(h, 'XDir', 'reverse');


end