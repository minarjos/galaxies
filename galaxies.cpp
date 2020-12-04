#include <vector>
#include <ncurses.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <set>
#include <chrono>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

// class for center of one galaxy
class center
{
public:
    // coordinates and index of center
    size_t x, y, n;
    // we keep centers in set, so they have to be comparable
    bool operator<(const center &other) const { return std::make_pair(x, y) < std::make_pair(other.x, other.y); };
    center(size_t x, size_t y, size_t n) : x(x), y(y), n(n) {}
    // reflect given square by this center, it is quite complicated because squares and centers use different coordinates
    std::pair<size_t, size_t> reflect(const std::pair<size_t, size_t> &pos) const { return {(2 * x - (2 * pos.first + 1)) / 2, (2 * y - (2 * pos.second + 1)) / 2}; }
};

// class for potential solution
class solution
{
    const int dx[4] = {-1, 0, 0, 1};
    const int dy[4] = {0, -1, 1, 0};
    size_t r, c;
    // we use dfs to check whether is given galaxy connected
    void dfs(int x, int y, size_t g, std::vector<std::vector<bool>> &visited)
    {
        if (x < 0 || x >= (int)r || y < 0 || y >= (int)c)
            return;
        if (visited[x][y] || galaxies[x][y] != g)
            return;
        visited[x][y] = true;
        for (size_t i = 0; i < 4; i++)
            dfs(x + dx[i], y + dy[i], g, visited);
    }

public:
    // galaxies[x][y] gives us index of galaxy in which is the square with coordinates x, y
    std::vector<std::vector<size_t>> galaxies;
    solution(size_t r, size_t c) : r(r), c(c) { galaxies = std::vector<std::vector<size_t>>(r, std::vector<size_t>(c, (size_t)0)); }
    // checks whether the solution is valid - every square is part of some galaxy and galaxies are connected
    bool is_valid()
    {
        std::vector<std::vector<bool>> visited = std::vector<std::vector<bool>>(r, std::vector<bool>(c, false));
        std::set<size_t> tested;
        // every square is in some galaxy
        for (auto &&row : galaxies)
            for (auto &&square : row)
                if (!square)
                    return false;
        // every galaxy is connected
        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < c; j++)
            {
                if (!visited[i][j] && tested.count(galaxies[i][j]))
                    return false;
                else if (!visited[i][j])
                {
                    dfs(i, j, galaxies[i][j], visited);
                    tested.insert(galaxies[i][j]);
                }
            }
        }
        return true;
    }
};

// stores possible next step of solution
struct possibility
{
    // index of galaxy
    size_t index;
    // distance from center
    size_t distance;
    // reflection of the square by the center
    std::pair<size_t, size_t> reflection;

    bool operator<(const possibility &other) const
    {
        return distance < other.distance;
    }
};

class board
{
    const int dx[4] = {-1, 0, 0, 1};
    const int dy[4] = {0, -1, 1, 0};
    //  returns true iff given square is in the board
    bool is_in(std::pair<size_t, size_t> pos) const
    {
        return pos.first >= 0 && pos.second >= 0 && pos.first < r && pos.second < c;
    }
    // we use dfs to check for reachable squares from given center
    void dfs(int x, int y, size_t g, const solution &s, std::vector<std::vector<bool>> &visited, std::vector<std::vector<std::set<size_t>>> &reachable) const
    {
        if (x < 0 || x >= (int)r || y < 0 || y >= (int)c)
            return;
        if (visited[x][y] || (s.galaxies[x][y] && s.galaxies[x][y] != g))
            return;
        reachable[x][y].insert(g);
        visited[x][y] = true;
        for (size_t i = 0; i < 4; i++)
            dfs(x + dx[i], y + dy[i], g, s, visited, reachable);
    }
    // for every square returns set of reachable centers
    std::vector<std::vector<std::set<size_t>>> get_reachable(const solution &s) const
    {
        std::vector<std::vector<std::set<size_t>>> reachable = std::vector<std::vector<std::set<size_t>>>(r, std::vector<std::set<size_t>>(c, std::set<size_t>()));
        for (auto &&cent : centers)
        {
            std::vector<std::vector<bool>> visited = std::vector<std::vector<bool>>(r, std::vector<bool>(c, false));
            dfs(cent.x / 2, cent.y / 2, cent.n, s, visited, reachable);
        }
        return std::move(reachable);
    }
    // returns all possible next steps using this square
    std::vector<possibility> get_possibilities(int i, int j, const solution &s, std::vector<std::vector<std::set<size_t>>> &reachable) const
    {
        std::vector<possibility> possibilites;
        // tries to reflect around every center
        for (auto &&cent : centers)
        {
            std::pair<size_t, size_t> reflected = cent.reflect({i, j});
            // reflection has to be in board and both the reflection and the initial square have to be reachable from given center
            if (is_in(reflected) && !s.galaxies[reflected.first][reflected.second] && reachable[i][j].count(cent.n) && reachable[reflected.first][reflected.second].count(cent.n))
            {
                possibility pos;
                size_t dist = abs(2 * (int)i + 1 - (int)cent.x) + abs(2 * (int)j + 1 - (int)cent.y);
                possibilites.push_back({cent.n, dist, reflected});
            }
        }
        return std::move(possibilites);
    }

public:
    size_t r, c, n;
    std::set<center> centers;
    // adds or removes given center (in case it was already in the set)
    void add_center(size_t x, size_t y)
    {
        center new_center = center(x, y, n);
        if (centers.count(new_center) == 0)
            centers.insert(center(x, y, n++));
        else
            centers.erase(new_center);
    }
    board(size_t r, size_t c) : r(r), c(c), n(1) {}

    // recursive function, tries to extend given solution
    solution solve(solution &s) const
    {
        // if the given solution is already invalid, doesnt vaste time and ends
        if (s.is_valid())
            return s;
        // the square with least possibilities
        std::pair<size_t, size_t> best_square;
        std::vector<possibility> best_possibilites;
        std::vector<std::vector<std::set<size_t>>> reachable = get_reachable(s);
        // we try every square which is not part of any galaxy
        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < c; j++)
            {
                if (s.galaxies[i][j])
                    continue;
                std::vector<possibility> possibilites = get_possibilities(i, j, s, reachable);
                // if there are no possibilities for some square, there is no solution
                if (possibilites.empty())
                    return s;
                // replaces the best square if suitable
                if (best_possibilites.empty() || possibilites.size() < best_possibilites.size())
                {
                    best_square = {i, j};
                    best_possibilites = possibilites;
                }
            }
        }
        // sort the possibilities by distance - closer distance is better (but this heuristic does actually make no big difference)
        std::sort(best_possibilites.begin(), best_possibilites.end());
        for (size_t i = 0; i < best_possibilites.size(); i++)
        {
            // set the square and its reflection
            s.galaxies[best_square.first][best_square.second] = best_possibilites[i].index;
            s.galaxies[best_possibilites[i].reflection.first][best_possibilites[i].reflection.second] = best_possibilites[i].index;
            solution new_solution = solve(s);
            // if we got valid solution, return it
            if (new_solution.is_valid())
                return new_solution;
            // undo and try next possibility
            s.galaxies[best_square.first][best_square.second] = 0;
            s.galaxies[best_possibilites[i].reflection.first][best_possibilites[i].reflection.second] = 0;
        }
        return s;
    }
    // tries to find a solution
    solution solve() const
    {
        // starts with empty solution
        solution s = solution(r, c);
        // starts with the obvious squares, reurns in case of contradiction
        for (auto &&cent : centers)
        {
            if (s.galaxies[cent.x / 2][cent.y / 2])
                return s;
            s.galaxies[cent.x / 2][cent.y / 2] = cent.n;
            if (cent.x % 2 == 0)
            {
                if (s.galaxies[cent.x / 2 - 1][cent.y / 2])
                    return s;
                s.galaxies[cent.x / 2 - 1][cent.y / 2] = cent.n;
            }
            if (cent.y % 2 == 0)
            {
                if (s.galaxies[cent.x / 2][cent.y / 2 - 1])
                    return s;
                s.galaxies[cent.x / 2][cent.y / 2 - 1] = cent.n;
            }
            if (cent.x % 2 == 0 && cent.y % 2 == 0)
            {
                if (s.galaxies[cent.x / 2 - 1][cent.y / 2 - 1])
                    return s;
                s.galaxies[cent.x / 2 - 1][cent.y / 2 - 1] = cent.n;
            }
        }
        return solve(s);
    }
    // first kind of row for string representation
    std::string row1() const
    {
        std::string s = "+";
        for (size_t i = 0; i < c; i++)
            s += " - +";
        return s;
    }
    // second kind of row for string representation
    std::string row2() const
    {
        std::string s = "|";
        for (size_t i = 0; i < c; i++)
            s += "   |";
        return s;
    }
    // strin representation of empty board
    std::vector<std::string> empty() const
    {
        std::string r1 = row1();
        std::string r2 = row2();
        std::vector<std::string> v;
        v.push_back(r1);
        for (size_t i = 0; i < r; i++)
        {
            v.push_back(r2);
            v.push_back(r1);
        }
        return v;
    }
    // adds centers to empty board
    std::vector<std::string> to_string() const
    {
        std::vector<std::string> rows = empty();
        for (auto &&c : centers)
            rows[c.x][2 * c.y] = 'o';
        return rows;
    }
    // creates string representation of solution / 
    std::vector<std::string> to_solution_string(const solution &s) const
    {
        std::vector<std::string> rows = empty();
        for (size_t i = 0; i < r; i++)
            for (size_t j = 0; j < c - 1; j++)
                if (s.galaxies[i][j] == s.galaxies[i][j + 1])
                    rows[2 * i + 1][4 * j + 4] = ' ';
        for (size_t i = 0; i < c; i++)
            for (size_t j = 0; j < r - 1; j++)
                if (s.galaxies[j][i] == s.galaxies[j + 1][i])
                    rows[2 * j + 2][4 * i + 2] = ' ';
        for (size_t i = 0; i < r - 1; i++)
            for (size_t j = 0; j < c - 1; j++)
                if (s.galaxies[i][j] == s.galaxies[i][j + 1] && s.galaxies[i + 1][j] == s.galaxies[i + 1][j + 1] && s.galaxies[i][j] == s.galaxies[i + 1][j])
                    rows[2 * i + 2][4 * j + 4] = ' ';
        for (auto &&c : centers)
            rows[c.x][2 * c.y] = 'o';
        return rows;
    }
    // deletes bad centers after resize
    void update_centers()
    {
        std::set<center> new_centers;
        for (auto &&cent : centers)
            if (cent.x < 2 * r && cent.y < 2 * c)
                new_centers.insert(cent);
        centers = new_centers;
    }
};

std::ostream &operator<<(std::ostream &os, const board &b)
{
    os << b.r << " " << b.c << std::endl;
    std::vector<std::string> rows = b.to_string();
    for (auto &&row : rows)
        os << row << std::endl;
    return os;
}

std::ostream &operator<<(std::ostream &os, const std::pair<board, solution> &p)
{
    std::vector<std::string> rows = p.first.to_solution_string(p.second);
    for (auto &&row : rows)
        os << row << std::endl;
    return os;
}

std::istream &operator>>(std::istream &is, board &b)
{
    b.centers.clear();
    if (!(is >> b.r >> b.c))
    {
        b.r = b.c = 7;
        return is;
    }
    std::string row;
    std::getline(is, row);
    for (size_t i = 0; i < 2 * b.r + 1; i++)
    {
        if (!(std::getline(is, row)) || row.size() < 4 * b.c + 1)
            return is;
        for (size_t j = 0; j < 4 * b.c + 1; j += 2)
            if (row[j] == 'o')
                b.add_center(i, j / 2);
    }
    return is;
}

// class representing the user interface using ncurses
class ui
{
    // position of cursor
    int x, y;
    // initializes ncurses
    void init_ncurses() const
    {
        initscr();
        keypad(stdscr, TRUE);
        start_color();
        init_color(COLOR_WHITE, 700, 700, 700);
        init_color(COLOR_GREEN, 0, 700, 0);
        init_color(COLOR_RED, 1000, 0, 0);
        init_pair(1, COLOR_GREEN, COLOR_BLACK);
        init_pair(2, COLOR_RED, COLOR_BLACK);
    }
    // prints using ncurses - adds collors
    void print(const std::vector<std::string> &rows) const
    {
        clear();
        for (auto &&row : rows)
        {
            for (auto &&square : row)
            {
                if (square == 'o')
                    attron(COLOR_PAIR(1));
                else if (square == 'x')
                    attron(COLOR_PAIR(2));
                printw("%c", square);
                if (square == 'o')
                    attroff(COLOR_PAIR(1));
                else if (square == 'x')
                    attroff(COLOR_PAIR(2));
            }
            printw("\n");
        }
        refresh();
    }
    // solves the board and prints solution
    void solve(const board &b) const
    {
        auto start = std::chrono::steady_clock::now();
        solution s = b.solve();
        auto end = std::chrono::steady_clock::now();
        print_solution(b, s);
        if (!s.is_valid())
            printw("No solution\n");
        else
            printw("Solved in %d ms\n", (std::chrono::duration_cast<std::chrono::milliseconds>(end - start)).count());
        refresh();
    }
    // loads board from file
    void load(board &b)
    {
        x = y = 1;
        char file_name[80];
        print_board(b);
        printw("File name: ");
        getstr(file_name);
        std::ifstream ifs(file_name);
        ifs >> b;
        print_board(b);
        printw("Loaded from %s\n", file_name);
        ifs.close();
    }
    // saves board to file
    void save(const board &b) const
    {
        char file_name[80];
        print_board(b);
        printw("File name: ");
        getstr(file_name);
        std::ofstream of(file_name);
        of << b << std::endl;
        print_board(b);
        printw("Saved to %s\n", file_name);
        of.close();
    }

public:
    void print_board(const board &b) const
    {
        std::vector<std::string> rows = b.to_string();
        rows[x][2 * y] = 'x';
        print(rows);
    }
    void print_solution(const board &b, const solution &s) const
    {
        std::vector<std::string> rows = b.to_solution_string(s);
        print(rows);
    }
    ui() : x(1), y(1)
    {
        init_ncurses();
    }
    ~ui()
    {
        endwin();
    }
    bool get_action(board &b)
    {
        int action = getch();
        // arrows move cursor
        if (action == KEY_LEFT)
            y = std::max(y - 1, 1);
        else if (action == KEY_RIGHT)
            y = std::min(y + 1, (int)(2 * b.c) - 1);
        else if (action == KEY_UP)
            x = std::max(x - 1, 1);
        else if (action == KEY_DOWN)
            x = std::min(x + 1, (int)(2 * b.r) - 1);
        // smace toggles center
        else if (action == ' ')
            b.add_center(x, y);
        // r,R and c,C for resize
        else if (action == 'r')
        {
            b.r = std::max((int)b.r - 1, 1);
            x = std::min(x, (int)(2 * b.r) - 1);
            b.update_centers();
        }
        else if (action == 'R')
            b.r++;
        else if (action == 'c')
        {
            b.c = std::max((int)b.c - 1, 1);
            y = std::min(y, (int)(2 * b.c) - 1);
            b.update_centers();
        }
        else if (action == 'C')
            b.c++;
        // l to load from file
        else if (action == 'l')
        {
            load(b);
            return true;
        }
        // s to save in file
        else if (action == 's')
        {
            save(b);
            return true;
        }
        // escape ends 
        else if (action == 27)
            return false;
        // backspace to delete all centers
        else if (action == 263)
            b.centers.clear();
        // enter to solve
        else if (action == 10)
        {
            solve(b);
            return true;
        }
        print_board(b);
        return true;
    }
};

int main(int argc, char **argv)
{
    bool interactive = false;

    int opt;
    board b = board(7, 7);

    // -i starts editor
    while ((opt = getopt(argc, argv, "i")) != -1)
    {
        switch (opt)
        {
        case 'i':
            interactive = true;
            break;
        default:
            std::cerr << "Usage:" << argv[0] << "[-i]\n"
                      << std::endl;
            return 1;
        }
    }
    // without -i, reads board from stdin and prints it to stdout
    if (!interactive)
    {
        std::cin >> b;
        solution s = b.solve();
        if (s.is_valid())
            std::cout << std::make_pair(b, s);
        else
            std::cout << "No solution" << std::endl;
    }
    else
    {
        ui u;
        u.print_board(b);
        while (u.get_action(b))
            ;
    }
    return 0;
}