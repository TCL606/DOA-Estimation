function peak_idxs = FindPeak(arr, detect_len, find_break)
    if (nargin < 3)
        find_break = 0;
    end
    peak_idxs = [];
    for j = detect_len + 1: 1: length(arr) - detect_len
        if sum(arr(j) > arr(j - detect_len: j - 1)) == detect_len && ...
                sum(arr(j) > arr(j + 1: j + detect_len)) == detect_len
            peak_idxs = [peak_idxs, j];
            if find_break == 1
                break
            end
        end
    end
end

