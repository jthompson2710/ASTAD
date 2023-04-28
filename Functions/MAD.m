function MatZero = MAD(X)
    Mat = X;
    sz = size(Mat);
    MatZero = zeros(sz);
    M = max(Mat);
    if M == 2
        Mat(Mat > 0) = 2
    end
end