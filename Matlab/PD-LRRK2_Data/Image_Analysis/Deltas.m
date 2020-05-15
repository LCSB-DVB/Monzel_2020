function [DeltaVector] = Deltas(VectorSorted)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    steps = length(VectorSorted)-1;

    for s = 1:steps
        DeltaVector(s) = abs(VectorSorted(s) - VectorSorted(s+1));
    end

end

