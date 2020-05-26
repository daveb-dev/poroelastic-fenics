



def normalize_solution(u):
    "Normalize u: return u divided by max(u)"

    u_array = u.vector()[:]

    u_max = np.max(np.abs(u_array))

    u_array /= u_max

    u.vector()[:] = u_array

    return u
