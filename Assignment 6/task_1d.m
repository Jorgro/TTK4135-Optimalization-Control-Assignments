Q = 1/2*[
    2 0
    0 2
];
R = 2/2;
A = [
    1 0.5
    0 1
];
B = [
    0.125
    0.5
];
[K,S,e] = dlqr(A,B,Q,R,[]); 
disp(K);
disp(S);
disp(e);