export unpack_dataframe

"""
    unpack_dataframe(df)
Take a dataframe `df` produced by the `produce_or_load` function
and unpack its savename parameters into distinct columns.
"""
function unpack_dataframe(df)
    # collect all column names in df
    cols = names(df) .|> Symbol
    # remove path from the column list
    filter!(c -> c ≠ :path, cols)
    # unpack the parameters from path into a DataFrame
    pars = DataFrame(map(p -> parse_savename(p)[2], df.path))
    # join parameters and data into a single DataFrame
    [pars df[:,cols]]
end
