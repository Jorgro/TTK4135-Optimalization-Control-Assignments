Q = 1/2 * [
    4 0
    0 4
];

A = [
    1 0.1
    -0.1 1-0.1
];

R = 1/2 * 1;

B = [
    0
    0.1
];

[K, S, e] = dlqr(A, B, Q, R, []);
disp(K);
disp(S);
disp(e);
