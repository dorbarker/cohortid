const program = """
#!/usr/bin/env julia

import CohortIdentification
CohortIdentification.main(ARGS)
"""
const default_dir = joinpath(homedir(), ".local/bin")

function get_user_path()::String
    # User selected path is currently to used

    prompt = "Installation directory (should be in \$PATH): [$default]"

    println(prompt)

    response = readline()

    install_path = if response == ""
        default_dir
    else
        response
    end

    joinpath(install_path, "cohortid")

end

function write_program(install_path::String)

    write(install_path, program)

    chmod(install_path, 0o775)

    println("Installed cohortid to $install_path")
end

function force_xdg_path()

    println("Ensuring $default_dir exists...")
    mkpath(default_dir)

    cohortid_path = joinpath(default_dir, "cohortid")
    write_program(cohortid_path)

end

function force_userprofile()

    cohortid_path = joinpath(ENV["localappdata"], "cohortid")
    write_program(cohortid_path)
end


function install()

    if Sys.iswindows()

        force_userprofile()

    else

        force_xdg_path()
    end
end

install()
