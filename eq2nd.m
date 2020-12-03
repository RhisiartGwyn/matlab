% es_eq2grado.m
% Script per il calcolo delle radici dell’ equazione di
% secondo grado ax^2+bx+c=0
a = input('inserisci il primo coefficiente a = ');
b = input('inserisci il secondo coefficiente b = ');
c = input('inserisci il terzo coefficiente c = ');

if a == 0
    disp('L''equazione è di primo grado')
    
    if b ~= 0
        x = -c/b;
        fprintf('La soluzione è %6.5f\n', x)
        
    else 
        fprintf('Equazione impossibile\n')
        
    end
        
 else
     D = b^2 - 4*a*c;

    if D < 0
        disp('Le radici non sono reali')

    elseif D == 0
        x1 = -b/(2*a);
        fprintf('La radice è unica e vale x=%6.5f \n', x1);

    else      
        x1 = (-b-sqrt(D))/(2*a);
        x2 = (-b+sqrt(D))/(2*a);
        fprintf('Le radici sono x1=%6.5f, x2=%6.5f\n',x1,x2);

     end    
end