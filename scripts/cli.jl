module ExampleCLI

export options, int_values, float_values

function options(args, defaults)
    iseven(length(args)) || throw(ArgumentError("arguments must be --key value pairs"))
    parsed = copy(defaults)
    for i in 1:2:length(args)
        flag = args[i]
        startswith(flag, "--") || throw(ArgumentError("expected --key, got $flag"))
        key = flag[3:end]
        haskey(defaults, key) || throw(ArgumentError("unknown option $flag"))
        parsed[key] = args[i + 1]
    end
    return parsed
end

int_values(value) = parse.(Int, split(value, ','))
float_values(value) = parse.(Float64, split(value, ','))

end
