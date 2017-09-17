/*
 * Copyright (c) 2017 Jacob Lifshay
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <limits>

constexpr long double pi = 3.141592653589793238462643383279502884197169399375105820974944592307L;

struct Vec2
{
    float x;
    float y;
    constexpr Vec2() noexcept : x(), y()
    {
    }
    constexpr explicit Vec2(float v) noexcept : x(v), y(v)
    {
    }
    constexpr Vec2(float x, float y) noexcept : x(x), y(y)
    {
    }
    constexpr Vec2 operator+() const noexcept
    {
        return *this;
    }
    constexpr Vec2 operator-() const noexcept
    {
        return Vec2(-x, -y);
    }
    friend constexpr Vec2 operator+(Vec2 a, Vec2 b) noexcept
    {
        return Vec2(a.x + b.x, a.y + b.y);
    }
    friend constexpr Vec2 operator-(Vec2 a, Vec2 b) noexcept
    {
        return Vec2(a.x - b.x, a.y - b.y);
    }
    friend constexpr Vec2 operator*(Vec2 a, Vec2 b) noexcept
    {
        return Vec2(a.x * b.x, a.y * b.y);
    }
    friend constexpr Vec2 operator*(float a, Vec2 b) noexcept
    {
        return Vec2(a * b.x, a * b.y);
    }
    friend constexpr Vec2 operator*(Vec2 a, float b) noexcept
    {
        return Vec2(a.x * b, a.y * b);
    }
    friend constexpr Vec2 operator/(Vec2 a, Vec2 b) noexcept
    {
        return Vec2(a.x / b.x, a.y / b.y);
    }
    friend constexpr Vec2 operator/(Vec2 a, float b) noexcept
    {
        return Vec2(a.x / b, a.y / b);
    }
    constexpr Vec2 &operator+=(Vec2 rt) noexcept
    {
        *this = *this + rt;
        return *this;
    }
    constexpr Vec2 &operator-=(Vec2 rt) noexcept
    {
        *this = *this - rt;
        return *this;
    }
    constexpr Vec2 &operator*=(Vec2 rt) noexcept
    {
        *this = *this * rt;
        return *this;
    }
    constexpr Vec2 &operator*=(float rt) noexcept
    {
        *this = *this * rt;
        return *this;
    }
    constexpr Vec2 &operator/=(Vec2 rt) noexcept
    {
        *this = *this / rt;
        return *this;
    }
    constexpr Vec2 &operator/=(float rt) noexcept
    {
        *this = *this / rt;
        return *this;
    }
    friend constexpr bool operator==(Vec2 a, Vec2 b) noexcept
    {
        return a.x == b.x && a.y == b.y;
    }
    friend constexpr bool operator!=(Vec2 a, Vec2 b) noexcept
    {
        return !(a == b);
    }
};

constexpr float dot(Vec2 a, Vec2 b) noexcept
{
    auto p = a * b;
    return p.x + p.y;
}

constexpr float abs_squared(Vec2 v) noexcept
{
    return dot(v, v);
}

inline float abs(Vec2 v) noexcept
{
    return std::sqrt(abs_squared(v));
}

inline Vec2 normalize(Vec2 v) noexcept
{
    float abs_squared_value = abs_squared(v);
    if(abs_squared_value == 0)
        return Vec2();
    return v / std::sqrt(abs_squared_value);
}

constexpr Vec2 interpolate(float t, Vec2 v0, Vec2 v1) noexcept
{
    return (1 - t) * v0 + t * v1;
}

constexpr float interpolate(float t, float v0, float v1) noexcept
{
    return (1 - t) * v0 + t * v1;
}

float angle_between(Vec2 a, Vec2 b) noexcept
{
    float denominator = std::sqrt(abs_squared(a) * abs_squared(b));
    if(denominator == 0) // a or b are zero
        return 0; // actual angle is not defined; pick zero
    float numerator = dot(a, b);
    float acos_arg = numerator / denominator;
    if(acos_arg < -1) // fix rounding errors
        acos_arg = -1;
    if(acos_arg > 1) // fix rounding errors
        acos_arg = 1;
    float positive_angle = std::acos(acos_arg);
    if(a.x * b.y < a.y * b.x)
        return -positive_angle;
    return positive_angle;
}

struct Cubic_curve
{
    Vec2 p0;
    Vec2 p1;
    Vec2 p2;
    Vec2 p3;
    constexpr Cubic_curve() noexcept : p0(), p1(), p2(), p3()
    {
    }
    constexpr explicit Cubic_curve(Vec2 v) noexcept : p0(v), p1(v), p2(v), p3(v)
    {
    }
    constexpr Cubic_curve(Vec2 p0, Vec2 p3) noexcept
        : p0(p0),
          p1(static_cast<float>(2 / 3.0) * p0 + static_cast<float>(1 / 3.0) * p3),
          p2(static_cast<float>(1 / 3.0) * p0 + static_cast<float>(2 / 3.0) * p3),
          p3(p3)
    {
    }
    constexpr Cubic_curve(Vec2 p0, Vec2 pc, Vec2 p3) noexcept
        : p0(p0),
          p1(static_cast<float>(1 / 3.0) * p0 + static_cast<float>(2 / 3.0) * pc),
          p2(static_cast<float>(2 / 3.0) * pc + static_cast<float>(1 / 3.0) * p3),
          p3(p3)
    {
    }
    constexpr Cubic_curve(Vec2 p0, Vec2 p1, Vec2 p2, Vec2 p3) noexcept : p0(p0),
                                                                         p1(p1),
                                                                         p2(p2),
                                                                         p3(p3)
    {
    }
    constexpr Vec2 evaluate(float t) const noexcept
    {
        float nt = 1 - t;
        return nt * nt * nt * p0 + nt * nt * t * 3.0f * p1 + nt * t * t * 3.0f * p2
               + t * t * t * p3;
    }
    constexpr Vec2 evaluate_first_derivative(float t) const noexcept
    {
        return 3.0f * ((p1 - p0)
                       + t * (2.0f * (p0 - p1 + (p2 - p1)) + t * (p3 - p0 + 3.0f * (p1 - p2))));
    }
    constexpr Vec2 evaluate_second_derivative(float t) const noexcept
    {
        return 6.0f * ((p2 - p1 + (p0 - p1)) + t * (p3 - p0 + 3.0f * (p1 - p2)));
    }
    constexpr Vec2 evaluate_third_derivative() const noexcept
    {
        return 6.0f * (p3 - p0 + 3.0f * (p1 - p2));
    }
    constexpr Cubic_curve reversed() const noexcept
    {
        return Cubic_curve(p3, p2, p1, p0);
    }
    constexpr Cubic_curve first_part() const noexcept
    {
        return Cubic_curve(p0,
                           (p1 + p0) * 0.5f,
                           (p2 + 2.0f * p1 + p0) * 0.25f,
                           (p3 + p0 + 3.0f * (p1 + p2)) * 0.125f);
    }
    constexpr Cubic_curve last_part() const noexcept
    {
        return reversed().first_part().reversed();
    }
    friend constexpr bool operator==(const Cubic_curve &a, const Cubic_curve &b) noexcept
    {
        return a.p0 == b.p0 && a.p1 == b.p1 && a.p2 == b.p2 && a.p3 == b.p3;
    }
    friend constexpr bool operator!=(const Cubic_curve &a, const Cubic_curve &b) noexcept
    {
        return !(a == b);
    }
};

struct Line
{
    Vec2 p0;
    Vec2 p1;
    constexpr Line() noexcept : p0(), p1()
    {
    }
    constexpr explicit Line(Vec2 v) noexcept : p0(v), p1(v)
    {
    }
    constexpr Line(Vec2 p0, Vec2 p1) noexcept : p0(p0), p1(p1)
    {
    }
    constexpr Vec2 evaluate(float t) const noexcept
    {
        float nt = 1 - t;
        return nt * p0 + t * p1;
    }
    constexpr Line reversed() const noexcept
    {
        return Line(p1, p0);
    }
    constexpr Line first_part() const noexcept
    {
        return Line(p0, (p1 + p0) * 0.5f);
    }
    constexpr Line last_part() const noexcept
    {
        return reversed().first_part().reversed();
    }
    friend constexpr bool operator==(const Line &a, const Line &b) noexcept
    {
        return a.p0 == b.p0 && a.p1 == b.p1;
    }
    friend constexpr bool operator!=(const Line &a, const Line &b) noexcept
    {
        return !(a == b);
    }
};

inline bool is_point_near_line(float tolerance, Vec2 p0, Vec2 p1, Vec2 p) noexcept
{
    float a = p0.y - p1.y;
    float b = p1.x - p0.x;
    float c = p0.x * p1.y - p1.x * p0.y;
    float adjusted_tolerance =
        tolerance * ((p1.x - p0.x) * p.y + (p.x - p1.x) * p0.y + (p0.x - p.x) * p1.y);
    return std::fabs(p.x * a + p.y * b + c) < std::fabs(adjusted_tolerance);
}

constexpr float convert_to_lines_default_tolerance = 1e-4;

inline std::vector<Line> convert_to_lines(Cubic_curve curve,
                                          std::vector<Line> prepended_lines = {},
                                          float tolerance = convert_to_lines_default_tolerance)
{
    bool need_subdivision{};
    if(abs_squared(curve.p0 - curve.p1) < tolerance * tolerance
       && abs_squared(curve.p0 - curve.p2) < tolerance * tolerance
       && abs_squared(curve.p0 - curve.p3) < tolerance * tolerance)
        need_subdivision = false;
    else
    {
        auto ideal_line = Cubic_curve(curve.p0, curve.p3);
        if(abs_squared(curve.p1 - ideal_line.p1) < tolerance * tolerance
           && abs_squared(curve.p2 - ideal_line.p2) < tolerance * tolerance)
            need_subdivision = false;
        else
            need_subdivision = true;
    }
    if(need_subdivision)
        return convert_to_lines(
            curve.last_part(),
            convert_to_lines(curve.first_part(), std::move(prepended_lines), tolerance),
            tolerance);
    prepended_lines.push_back(Line(curve.p0, curve.p3));
    return prepended_lines;
}

inline std::vector<Line> convert_to_lines(Cubic_curve curve, float tolerance)
{
    return convert_to_lines(curve, {}, tolerance);
}

inline std::vector<Line> convert_to_lines(const Cubic_curve *curves,
                                          std::size_t curve_count,
                                          std::vector<Line> prepended_lines = {},
                                          float tolerance = convert_to_lines_default_tolerance)
{
    for(std::size_t i = 0; i < curve_count; i++)
        prepended_lines = convert_to_lines(curves[i], std::move(prepended_lines), tolerance);
    return prepended_lines;
}

inline std::vector<Line> convert_to_lines(const Cubic_curve *curves,
                                          std::size_t curve_count,
                                          float tolerance)
{
    return convert_to_lines(curves, curve_count, {}, tolerance);
}

struct Svg_parse_error : public std::runtime_error
{
    std::size_t line_number;
    std::size_t line_byte;
    static std::string make_message(std::size_t line_number,
                                    std::size_t line_byte,
                                    const std::string &source_name,
                                    const std::string &error_message)
    {
        std::ostringstream ss;
        ss << source_name << ":" << line_number << ":" << line_byte << ": error: " << error_message;
        return ss.str();
    }
    Svg_parse_error(std::size_t line_number,
                    std::size_t line_byte,
                    const std::string &source_name,
                    const std::string &error_message)
        : std::runtime_error(make_message(line_number, line_byte, source_name, error_message)),
          line_number(line_number),
          line_byte(line_byte)
    {
    }
};

class Svg_path_interpreter
{
private:
    std::vector<Cubic_curve> curves;
    Vec2 close_path_target;
    bool has_close_path_target = false;
    Vec2 current_position = Vec2(0);
    Vec2 last_cubic_control_point = Vec2(0);
    Vec2 last_quadratic_control_point = Vec2(0);

private:
    void set_current_position(Vec2 position) noexcept
    {
        current_position = position;
        last_cubic_control_point = position;
        last_quadratic_control_point = position;
    }
    void set_current_position_and_last_cubic_control_point(Vec2 position,
                                                           Vec2 control_point) noexcept
    {
        current_position = position;
        last_cubic_control_point = control_point;
        last_quadratic_control_point = position;
    }
    void set_current_position_and_last_quadratic_control_point(Vec2 position,
                                                               Vec2 control_point) noexcept
    {
        current_position = position;
        last_cubic_control_point = position;
        last_quadratic_control_point = control_point;
    }
    /** doesn't set the current position */
    void append_elliptical_arc(Vec2 start_position,
                               Vec2 end_position,
                               Vec2 radius,
                               float sin_x_axis_rotation,
                               float cos_x_axis_rotation,
                               Vec2 center_position,
                               float start_angle_in_radians,
                               float delta_angle_in_radians)
    {
        constexpr unsigned curves_per_circle = 12;
        constexpr float curves_per_radian = curves_per_circle / (2 * pi);
        unsigned curve_count =
            static_cast<unsigned>(std::ceil(std::fabs(delta_angle_in_radians) * curves_per_radian));
        if(curve_count > curves_per_circle) // prevent infinite loops due to NaN or Infinity
            curve_count = curves_per_circle;
        else if(curve_count == 0) // we need at least 1 curve
            curve_count = 1;
        float delta_radians_per_curve = delta_angle_in_radians / static_cast<float>(curve_count);
        auto get_position_on_arc = [&](Vec2 cossin_angle) noexcept
        {
            Vec2 untransformed_position = radius * cossin_angle;
            Vec2 untranslated_position(cos_x_axis_rotation * untransformed_position.x
                                           - sin_x_axis_rotation * untransformed_position.y,
                                       cos_x_axis_rotation * untransformed_position.y
                                           + sin_x_axis_rotation * untransformed_position.x);
            return center_position + untranslated_position;
        };
        auto get_unit_circle_arc_first_control_point = [](Vec2 cossin_start_angle,
                                                          Vec2 cossin_end_angle) noexcept
        {
            Vec2 cossin_center_angle =
                normalize(cossin_start_angle
                          + cossin_end_angle); // get the position for the angle halfway in between
            return Vec2((cossin_start_angle.y * (4 * cossin_end_angle.x * cossin_end_angle.x
                                                 - 8 * cossin_end_angle.x * cossin_center_angle.x
                                                 + 4 * cossin_end_angle.y * cossin_end_angle.y
                                                 - 8 * cossin_end_angle.y * cossin_center_angle.y)
                         + 4 * cossin_start_angle.y * cossin_start_angle.y * cossin_end_angle.y
                         + 3 * cossin_start_angle.x * cossin_start_angle.x * cossin_end_angle.y
                         + cossin_start_angle.x * cossin_end_angle.x * cossin_start_angle.y)
                            / (3 * cossin_start_angle.x * cossin_end_angle.y
                               - 3 * cossin_end_angle.x * cossin_start_angle.y),
                        -(cossin_start_angle.x * (4 * cossin_end_angle.x * cossin_end_angle.x
                                                  - 8 * cossin_end_angle.x * cossin_center_angle.x
                                                  + cossin_start_angle.y * cossin_end_angle.y
                                                  + 4 * cossin_end_angle.y * cossin_end_angle.y
                                                  - 8 * cossin_end_angle.y * cossin_center_angle.y)
                          + 3 * cossin_end_angle.x * cossin_start_angle.y * cossin_start_angle.y
                          + 4 * cossin_start_angle.x * cossin_start_angle.x * cossin_end_angle.x)
                            / (3 * cossin_start_angle.x * cossin_end_angle.y
                               - 3 * cossin_end_angle.x * cossin_start_angle.y));
        };
        auto get_first_control_point_for_arc = [&](Vec2 cossin_start_angle,
                                                   Vec2 cossin_end_angle) noexcept
        {
            Vec2 untransformed_position = radius * get_unit_circle_arc_first_control_point(
                                                       cossin_start_angle, cossin_end_angle);
            Vec2 untranslated_position(cos_x_axis_rotation * untransformed_position.x
                                           - sin_x_axis_rotation * untransformed_position.y,
                                       cos_x_axis_rotation * untransformed_position.y
                                           + sin_x_axis_rotation * untransformed_position.x);
            return center_position + untranslated_position;
        };
        Vec2 curve_start_position{};
        // use the exact start position to prevent gaps/overlap from round-off errors
        Vec2 curve_end_position = start_position;
        float curve_start_angle{};
        float curve_end_angle = start_angle_in_radians;
        Vec2 cossin_curve_start_angle{};
        Vec2 cossin_curve_end_angle(std::cos(curve_end_angle), std::sin(curve_end_angle));
        for(unsigned curve_index = 0; curve_index < curve_count; curve_index++)
        {
            curve_start_angle = curve_end_angle;
            curve_end_angle += delta_radians_per_curve;
            cossin_curve_start_angle = cossin_curve_end_angle;
            cossin_curve_end_angle = Vec2(std::cos(curve_end_angle), std::sin(curve_end_angle));
            curve_start_position = curve_end_position;
            if(curve_index + 1 == curve_count)
            {
                // use the exact end position to prevent gaps/overlap from round-off errors
                curve_end_position = end_position;
            }
            else
            {
                curve_end_position = get_position_on_arc(cossin_curve_end_angle);
            }
            curves.push_back(Cubic_curve(
                curve_start_position,
                get_first_control_point_for_arc(cossin_curve_start_angle, cossin_curve_end_angle),
                get_first_control_point_for_arc(cossin_curve_end_angle, cossin_curve_start_angle),
                curve_end_position));
        }
    }

public:
    void move_to_absolute(Vec2 position) noexcept
    {
        if(!has_close_path_target)
            close_path_target = position;
        has_close_path_target = true;
        set_current_position(position);
    }
    void move_to_relative(Vec2 offset) noexcept
    {
        move_to_absolute(offset + current_position);
    }
    void close_path()
    {
        curves.push_back(Cubic_curve(current_position, close_path_target));
        has_close_path_target = false;
        set_current_position(close_path_target);
    }
    void line_to_absolute(Vec2 position)
    {
        curves.push_back(Cubic_curve(current_position, position));
        set_current_position(position);
    }
    void line_to_relative(Vec2 offset)
    {
        line_to_absolute(offset + current_position);
    }
    void horizontal_line_to_absolute(float x)
    {
        auto target = current_position;
        target.x = x;
        line_to_absolute(target);
    }
    void horizontal_line_to_relative(float x)
    {
        auto target = current_position;
        target.x += x;
        line_to_absolute(target);
    }
    void vertical_line_to_absolute(float y)
    {
        auto target = current_position;
        target.y = y;
        line_to_absolute(target);
    }
    void vertical_line_to_relative(float y)
    {
        auto target = current_position;
        target.y += y;
        line_to_absolute(target);
    }
    void cubic_curve_to_absolute(Vec2 p1, Vec2 p2, Vec2 p)
    {
        curves.push_back(Cubic_curve(current_position, p1, p2, p));
        set_current_position_and_last_cubic_control_point(p, p2);
    }
    void cubic_curve_to_relative(Vec2 p1, Vec2 p2, Vec2 p)
    {
        cubic_curve_to_absolute(p1 + current_position, p2 + current_position, p + current_position);
    }
    void cubic_shorthand_curve_to_absolute(Vec2 p2, Vec2 p)
    {
        cubic_curve_to_absolute(2 * current_position - last_cubic_control_point, p2, p);
    }
    void cubic_shorthand_curve_to_relative(Vec2 p2, Vec2 p)
    {
        cubic_shorthand_curve_to_absolute(p2 + current_position, p + current_position);
    }
    void quadratic_curve_to_absolute(Vec2 cp, Vec2 p)
    {
        curves.push_back(Cubic_curve(current_position, cp, p));
        set_current_position_and_last_quadratic_control_point(p, cp);
    }
    void quadratic_curve_to_relative(Vec2 cp, Vec2 p)
    {
        quadratic_curve_to_absolute(cp + current_position, p + current_position);
    }
    void quadratic_shorthand_curve_to_absolute(Vec2 p)
    {
        quadratic_curve_to_absolute(2 * current_position - last_quadratic_control_point, p);
    }
    void quadratic_shorthand_curve_to_relative(Vec2 p)
    {
        quadratic_shorthand_curve_to_absolute(p + current_position);
    }
    void elliptical_arc_to_absolute(
        Vec2 radius, float x_axis_rotation, bool large_arc_flag, bool sweep_flag, Vec2 end_position)
    {
        // algorithm from https://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes
        radius.x = std::fabs(radius.x);
        radius.y = std::fabs(radius.y);
        x_axis_rotation = std::fmod(x_axis_rotation, 360.0f);
        if(x_axis_rotation > 180)
            x_axis_rotation -= 360;
        if(x_axis_rotation <= -180)
            x_axis_rotation += 360;
        constexpr float degrees_to_radians = pi / 180;
        float x_axis_rotation_in_radians = degrees_to_radians * x_axis_rotation;
        float sin_x_axis_rotation = std::sin(x_axis_rotation_in_radians);
        float cos_x_axis_rotation = std::cos(x_axis_rotation_in_radians);
        auto start_position = current_position;
        auto translated_start_position = 0.5f * (start_position - end_position);
        Vec2 transformed_start_position(cos_x_axis_rotation * translated_start_position.x
                                            + sin_x_axis_rotation * translated_start_position.y,
                                        cos_x_axis_rotation * translated_start_position.y
                                            - sin_x_axis_rotation * translated_start_position.x);
        float lambda = abs_squared(transformed_start_position / radius);
        if(lambda > 1)
            radius *= std::sqrt(lambda);
        float radicand_denominator =
            abs_squared(Vec2(radius.y, radius.x) * transformed_start_position);
        float radicand_numerator =
            radius.x * radius.x * (radius.y * radius.y) - radicand_denominator;
        if(radicand_numerator == 0 || radicand_denominator == 0 || radius.x == 0 || radius.y == 0)
        {
            // either radius is zero or the arc has zero length
            line_to_absolute(end_position);
            return;
        }
        float radicand = radicand_numerator / radicand_denominator;
        if(radicand < 0) // can only happen because of round-off error
            radicand = 0;
        Vec2 transformed_unscaled_center_position(
            radius.x * transformed_start_position.y / radius.y,
            -radius.y * transformed_start_position.x / radius.x);
        Vec2 transformed_center_position =
            (large_arc_flag == sweep_flag ? -std::sqrt(radicand) : std::sqrt(radicand))
            * transformed_unscaled_center_position;
        Vec2 translated_center_position(cos_x_axis_rotation * transformed_center_position.x
                                            - sin_x_axis_rotation * transformed_center_position.y,
                                        cos_x_axis_rotation * transformed_center_position.y
                                            + sin_x_axis_rotation * transformed_center_position.x);
        Vec2 center_position = transformed_center_position + 0.5f * (start_position + end_position);
        auto angle_between_arg =
            (transformed_start_position - transformed_center_position) / radius;
        float start_angle_in_radians = angle_between(Vec2(1, 0), angle_between_arg);
        float delta_angle_in_radians =
            angle_between((transformed_start_position - transformed_center_position) / radius,
                          -(transformed_start_position + transformed_center_position) / radius);
        if(!sweep_flag && delta_angle_in_radians > 0)
            delta_angle_in_radians -= static_cast<float>(2 * pi);
        else if(sweep_flag && delta_angle_in_radians < 0)
            delta_angle_in_radians += static_cast<float>(2 * pi);
        append_elliptical_arc(start_position,
                              end_position,
                              radius,
                              sin_x_axis_rotation,
                              cos_x_axis_rotation,
                              center_position,
                              start_angle_in_radians,
                              delta_angle_in_radians);
        set_current_position(end_position);
    }
    void elliptical_arc_to_relative(
        Vec2 radius, float x_axis_rotation, bool large_arc_flag, bool sweep_flag, Vec2 end_position)
    {
        elliptical_arc_to_absolute(
            radius, x_axis_rotation, large_arc_flag, sweep_flag, end_position + current_position);
    }
    std::vector<Cubic_curve> finish() && noexcept
    {
        return std::move(curves);
    }
};

std::vector<Cubic_curve> parse_svg_path(const char *source,
                                        std::size_t source_size,
                                        const std::string &source_name)
{
    constexpr int eof = -1;
    struct Parser
    {
        int peek = ' ';
        std::size_t line_number = 1;
        std::size_t line_byte = 0;
        const char *source;
        std::size_t source_size;
        const std::string &source_name;
        void next() noexcept
        {
            if(source_size != 0)
            {
                int last_char = peek;
                peek = static_cast<unsigned char>(*source);
                source++;
                source_size--;
                if(last_char == '\r' && peek == '\n')
                {
                    line_byte++;
                }
                else if(last_char == '\n' || last_char == '\r')
                {
                    line_byte = 1;
                    line_number++;
                }
                else
                {
                    line_byte++;
                }
            }
            else
                peek = eof;
        };
        int get() noexcept
        {
            int retval = peek;
            next();
            return retval;
        };

        void parse_error(const std::string &error_message)
        {
            throw Svg_parse_error(line_number, line_byte, source_name, error_message);
        };

        Svg_path_interpreter path_interpreter;

        // parser based off of BNF grammar in https://www.w3.org/TR/SVG/paths.html#PathDataBNF
        bool is_wsp(int ch) noexcept
        {
            switch(ch)
            {
            case ' ':
            case '\t':
            case '\r':
            case '\n':
                return true;
            }
            return false;
        };
        bool is_digit(int ch) noexcept
        {
            if(ch >= '0' && ch <= '9')
                return true;
            return false;
        };
        bool is_sign(int ch) noexcept
        {
            if(ch == '+' || ch == '-')
                return true;
            return false;
        };
        bool is_exponent_start(int ch) noexcept
        {
            if(ch == 'e' || ch == 'E')
                return true;
            return false;
        };
        bool is_flag(int ch) noexcept
        {
            if(ch == '0' || ch == '1')
                return true;
            return false;
        };
        bool is_nonnegative_number_start(int ch) noexcept
        {
            if(ch == '.' || is_digit(ch))
                return true;
            return false;
        };
        bool is_number_start(int ch) noexcept
        {
            if(is_nonnegative_number_start(ch) || is_sign(ch))
                return true;
            return false;
        };
        bool is_comma(int ch) noexcept
        {
            if(ch == ',')
                return true;
            return false;
        };
        bool is_comma_wsp_start(int ch) noexcept
        {
            if(ch == ',' || is_wsp(ch))
                return true;
            return false;
        };
        bool is_coordinate_start(int ch) noexcept
        {
            return is_number_start(ch);
        };
        bool is_coordinate_pair_start(int ch) noexcept
        {
            return is_coordinate_start(ch);
        };
        bool is_elliptical_arc_argument_start(int ch) noexcept
        {
            return is_nonnegative_number_start(ch);
        };
        bool is_elliptical_arc_argument_sequence_start(int ch) noexcept
        {
            return is_elliptical_arc_argument_start(ch);
        };
        bool is_elliptical_arc_start(int ch) noexcept
        {
            if(ch == 'a' || ch == 'A')
                return true;
            return false;
        };
        bool is_smooth_quadratic_bezier_curveto_argument_sequence_start(int ch) noexcept
        {
            return is_coordinate_pair_start(ch);
        };
        bool is_smooth_quadratic_bezier_curveto_start(int ch) noexcept
        {
            if(ch == 't' || ch == 'T')
                return true;
            return false;
        };
        bool is_quadratic_bezier_curveto_argument_start(int ch) noexcept
        {
            return is_coordinate_pair_start(ch);
        };
        bool is_quadratic_bezier_curveto_argument_sequence_start(int ch) noexcept
        {
            return is_quadratic_bezier_curveto_argument_start(ch);
        };
        bool is_quadratic_bezier_curveto_start(int ch) noexcept
        {
            if(ch == 'q' || ch == 'Q')
                return true;
            return false;
        };
        bool is_smooth_curveto_argument_start(int ch) noexcept
        {
            return is_coordinate_pair_start(ch);
        };
        bool is_smooth_curveto_argument_sequence_start(int ch) noexcept
        {
            return is_smooth_curveto_argument_start(ch);
        };
        bool is_smooth_curveto_start(int ch) noexcept
        {
            if(ch == 's' || ch == 'S')
                return true;
            return false;
        };
        bool is_curveto_argument_start(int ch) noexcept
        {
            return is_coordinate_pair_start(ch);
        };
        bool is_curveto_argument_sequence_start(int ch) noexcept
        {
            return is_curveto_argument_start(ch);
        };
        bool is_curveto_start(int ch) noexcept
        {
            if(ch == 'c' || ch == 'C')
                return true;
            return false;
        };
        bool is_vertical_lineto_argument_sequence_start(int ch) noexcept
        {
            return is_coordinate_start(ch);
        };
        bool is_vertical_lineto_start(int ch) noexcept
        {
            if(ch == 'v' || ch == 'V')
                return true;
            return false;
        };
        bool is_horizontal_lineto_argument_sequence_start(int ch) noexcept
        {
            return is_coordinate_start(ch);
        };
        bool is_horizontal_lineto_start(int ch) noexcept
        {
            if(ch == 'h' || ch == 'H')
                return true;
            return false;
        };
        bool is_lineto_argument_sequence_start(int ch) noexcept
        {
            return is_coordinate_pair_start(ch);
        };
        bool is_lineto_start(int ch) noexcept
        {
            if(ch == 'l' || ch == 'L')
                return true;
            return false;
        };
        bool is_closepath_start(int ch) noexcept
        {
            if(ch == 'z' || ch == 'Z')
                return true;
            return false;
        };
        bool is_moveto_argument_sequence_start(int ch) noexcept
        {
            return is_coordinate_pair_start(ch);
        };
        bool is_moveto_start(int ch) noexcept
        {
            if(ch == 'm' || ch == 'M')
                return true;
            return false;
        };
        bool is_drawto_command_start(int ch) noexcept
        {
            if(is_closepath_start(ch))
                return true;
            if(is_lineto_start(ch))
                return true;
            if(is_horizontal_lineto_start(ch))
                return true;
            if(is_vertical_lineto_start(ch))
                return true;
            if(is_curveto_start(ch))
                return true;
            if(is_smooth_curveto_start(ch))
                return true;
            if(is_quadratic_bezier_curveto_start(ch))
                return true;
            if(is_smooth_quadratic_bezier_curveto_start(ch))
                return true;
            if(is_elliptical_arc_start(ch))
                return true;
            return false;
        };
        bool is_drawto_commands_start(int ch) noexcept
        {
            return is_drawto_command_start(ch);
        };
        bool is_moveto_drawto_command_group_start(int ch) noexcept
        {
            return is_moveto_start(ch);
        };
        bool is_moveto_drawto_command_groups_start(int ch) noexcept
        {
            return is_moveto_drawto_command_group_start(ch);
        };

        void parse_wsp() noexcept
        {
            assert(is_wsp(peek));
            next();
        };
        int parse_digit() noexcept
        {
            assert(is_digit(peek));
            return get() - '0';
        };
        struct Digit_sequence
        {
            std::size_t parsed_digit_count = 0; // number of digits parsed
            std::size_t kept_digit_count = 0; // number of digits in value skipping initial zeros
            std::size_t ignored_digit_count = 0; // number of digits after value
            std::uint64_t value = 0;
            void append_digit(int digit) noexcept
            {
                parsed_digit_count++;
                constexpr std::size_t max_kept_digit_count = 17;
                if(kept_digit_count < max_kept_digit_count && (value != 0 || digit != 0))
                {
                    value *= 10;
                    value += digit;
                    kept_digit_count++;
                }
                else if(kept_digit_count != 0)
                    ignored_digit_count++;
            }
        };
        Digit_sequence parse_digit_sequence(Digit_sequence digit_sequence) noexcept
        {
            assert(is_digit(peek));
            while(is_digit(peek))
                digit_sequence.append_digit(parse_digit());
            return digit_sequence;
        };
        bool parse_sign() noexcept
        {
            assert(is_sign(peek));
            return get() == '-';
        };
        std::int64_t parse_exponent()
        {
            assert(is_exponent_start(peek));
            next();
            bool sign = is_sign(peek) && parse_sign();
            if(!is_digit(peek))
                parse_error("missing digit in exponent");
            auto digit_sequence = parse_digit_sequence({});
            if(digit_sequence.ignored_digit_count != 0) // too many digits, just return a big number
                return sign ? -(1LL << 60) : 1LL << 60;
            auto retval = static_cast<std::int64_t>(digit_sequence.value);
            if(sign)
                retval = -retval;
            return retval;
        };
        double parse_nonnegative_number()
        {
            assert(is_nonnegative_number_start(peek));
            Digit_sequence mantissa;
            if(is_digit(peek))
                mantissa = parse_digit_sequence(mantissa);
            std::size_t digits_before_dot = mantissa.parsed_digit_count;
            if(peek == '.')
            {
                next();
                if(is_digit(peek))
                    mantissa = parse_digit_sequence(mantissa);
            }
            if(mantissa.parsed_digit_count == 0)
                parse_error("missing digit in number");
            std::ptrdiff_t digits_after_dot = mantissa.parsed_digit_count - digits_before_dot;
            std::ptrdiff_t digits_in_value_after_dot =
                digits_after_dot - static_cast<std::ptrdiff_t>(mantissa.ignored_digit_count);
            std::int64_t exponent = -digits_in_value_after_dot;
            if(is_exponent_start(peek))
                exponent += parse_exponent();
            double retval = mantissa.value;
            if(exponent < -100)
            {
                // avoid incorrectly rounding to zero for large negative exponents
                retval *= 1e-100;
                exponent += 100;
            }
            return retval * std::pow(10.0, exponent);
        };
        double parse_number()
        {
            assert(is_number_start(peek));
            bool sign = is_sign(peek) && parse_sign();
            if(!is_nonnegative_number_start(peek))
                parse_error("missing digit or '.' in number");
            double retval = parse_nonnegative_number();
            if(sign)
                return -retval;
            return retval;
        };
        bool parse_flag() noexcept
        {
            assert(is_flag(peek));
            return get() != '0';
        };
        void parse_comma_wsp() noexcept
        {
            assert(is_comma_wsp_start(peek));
            while(is_wsp(peek))
                parse_wsp();
            if(is_comma(peek))
                next();
            while(is_wsp(peek))
                parse_wsp();
        };
        double parse_coordinate()
        {
            assert(is_coordinate_start(peek));
            return parse_number();
        };
        Vec2 parse_coordinate_pair()
        {
            assert(is_coordinate_pair_start(peek));
            double x = parse_coordinate();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_coordinate_start(peek))
                parse_error("missing y coordinate");
            double y = parse_coordinate();
            return Vec2(x, y);
        };
        void parse_elliptical_arc_argument(bool is_relative)
        {
            assert(is_elliptical_arc_argument_start(peek));
            float radius_x = parse_nonnegative_number();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_nonnegative_number_start(peek))
                parse_error("elliptical arc missing radius y component");
            float radius_y = parse_nonnegative_number();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_number_start(peek))
                parse_error("elliptical arc missing x axis rotation");
            float x_axis_rotation = parse_number();
            if(!is_comma_wsp_start(peek))
                parse_error("missing comma or whitespace");
            parse_comma_wsp();
            if(!is_flag(peek))
                parse_error("missing large arc flag");
            bool large_arc_flag = parse_flag();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_flag(peek))
                parse_error("missing sweep flag");
            bool sweep_flag = parse_flag();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_coordinate_pair_start(peek))
                parse_error("elliptical arc missing end position");
            auto end_position = parse_coordinate_pair();
            if(is_relative)
                path_interpreter.elliptical_arc_to_relative(Vec2(radius_x, radius_y),
                                                            x_axis_rotation,
                                                            large_arc_flag,
                                                            sweep_flag,
                                                            end_position);
            else
                path_interpreter.elliptical_arc_to_absolute(Vec2(radius_x, radius_y),
                                                            x_axis_rotation,
                                                            large_arc_flag,
                                                            sweep_flag,
                                                            end_position);
        };
        void parse_elliptical_arc_argument_sequence(bool is_relative)
        {
            assert(is_elliptical_arc_argument_sequence_start(peek));
            while(true)
            {
                parse_elliptical_arc_argument(is_relative);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_elliptical_arc_argument_start(peek))
                    break;
            }
        };
        void parse_elliptical_arc()
        {
            assert(is_elliptical_arc_start(peek));
            bool is_relative = get() == 'a';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_elliptical_arc_argument_sequence_start(peek))
                parse_error("missing elliptical arc argument sequence");
            parse_elliptical_arc_argument_sequence(is_relative);
        };
        void parse_smooth_quadratic_bezier_curveto_argument_sequence(bool is_relative)
        {
            assert(is_smooth_quadratic_bezier_curveto_argument_sequence_start(peek));
            while(true)
            {
                auto position = parse_coordinate_pair();
                if(is_relative)
                    path_interpreter.quadratic_shorthand_curve_to_relative(position);
                else
                    path_interpreter.quadratic_shorthand_curve_to_absolute(position);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_coordinate_pair_start(peek))
                    break;
            }
        };
        void parse_smooth_quadratic_bezier_curveto()
        {
            assert(is_smooth_quadratic_bezier_curveto_start(peek));
            bool is_relative = get() == 't';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_smooth_quadratic_bezier_curveto_argument_sequence_start(peek))
                parse_error("missing smooth quadratic bezier curveto argument sequence");
            parse_smooth_quadratic_bezier_curveto_argument_sequence(is_relative);
        };
        void parse_quadratic_bezier_curveto_argument(bool is_relative)
        {
            assert(is_quadratic_bezier_curveto_argument_start(peek));
            auto control_point = parse_coordinate_pair();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_coordinate_pair_start(peek))
                parse_error("quadratic bezier curveto missing end position");
            auto end_position = parse_coordinate_pair();
            if(is_relative)
                path_interpreter.quadratic_curve_to_relative(control_point, end_position);
            else
                path_interpreter.quadratic_curve_to_absolute(control_point, end_position);
        };
        void parse_quadratic_bezier_curveto_argument_sequence(bool is_relative)
        {
            assert(is_quadratic_bezier_curveto_argument_sequence_start(peek));
            while(true)
            {
                parse_quadratic_bezier_curveto_argument(is_relative);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_quadratic_bezier_curveto_argument_start(peek))
                    break;
            }
        };
        void parse_quadratic_bezier_curveto()
        {
            assert(is_quadratic_bezier_curveto_start(peek));
            bool is_relative = get() == 'q';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_quadratic_bezier_curveto_argument_sequence_start(peek))
                parse_error("missing quadratic bezier curveto argument sequence");
            parse_quadratic_bezier_curveto_argument_sequence(is_relative);
        };
        void parse_smooth_curveto_argument(bool is_relative)
        {
            assert(is_smooth_curveto_argument_start(peek));
            auto control_point = parse_coordinate_pair();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_coordinate_pair_start(peek))
                parse_error("smooth curveto missing end position");
            auto end_position = parse_coordinate_pair();
            if(is_relative)
                path_interpreter.cubic_shorthand_curve_to_relative(control_point, end_position);
            else
                path_interpreter.cubic_shorthand_curve_to_absolute(control_point, end_position);
        };
        void parse_smooth_curveto_argument_sequence(bool is_relative)
        {
            assert(is_smooth_curveto_argument_sequence_start(peek));
            while(true)
            {
                parse_smooth_curveto_argument(is_relative);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_smooth_curveto_argument_start(peek))
                    break;
            }
        };
        void parse_smooth_curveto()
        {
            assert(is_smooth_curveto_start(peek));
            bool is_relative = get() == 's';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_smooth_curveto_argument_sequence_start(peek))
                parse_error("missing smooth curveto argument sequence");
            parse_smooth_curveto_argument_sequence(is_relative);
        };
        void parse_curveto_argument(bool is_relative)
        {
            assert(is_curveto_argument_start(peek));
            auto first_control_point = parse_coordinate_pair();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_coordinate_pair_start(peek))
                parse_error("curveto missing second control point");
            auto second_control_point = parse_coordinate_pair();
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(!is_coordinate_pair_start(peek))
                parse_error("curveto missing end position");
            auto end_position = parse_coordinate_pair();
            if(is_relative)
                path_interpreter.cubic_curve_to_relative(
                    first_control_point, second_control_point, end_position);
            else
                path_interpreter.cubic_curve_to_absolute(
                    first_control_point, second_control_point, end_position);
        };
        void parse_curveto_argument_sequence(bool is_relative)
        {
            assert(is_curveto_argument_sequence_start(peek));
            while(true)
            {
                parse_curveto_argument(is_relative);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_curveto_argument_start(peek))
                    break;
            }
        };
        void parse_curveto()
        {
            assert(is_curveto_start(peek));
            bool is_relative = get() == 'c';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_curveto_argument_sequence_start(peek))
                parse_error("missing curveto argument sequence");
            parse_curveto_argument_sequence(is_relative);
        };
        void parse_vertical_lineto_argument_sequence(bool is_relative)
        {
            assert(is_vertical_lineto_argument_sequence_start(peek));
            while(true)
            {
                auto position = parse_coordinate();
                if(is_relative)
                    path_interpreter.vertical_line_to_relative(position);
                else
                    path_interpreter.vertical_line_to_absolute(position);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_coordinate_start(peek))
                    break;
            }
        };
        void parse_vertical_lineto()
        {
            assert(is_vertical_lineto_start(peek));
            bool is_relative = get() == 'v';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_vertical_lineto_argument_sequence_start(peek))
                parse_error("missing vertical lineto argument sequence");
            parse_vertical_lineto_argument_sequence(is_relative);
        };
        void parse_horizontal_lineto_argument_sequence(bool is_relative)
        {
            assert(is_horizontal_lineto_argument_sequence_start(peek));
            while(true)
            {
                auto position = parse_coordinate();
                if(is_relative)
                    path_interpreter.horizontal_line_to_relative(position);
                else
                    path_interpreter.horizontal_line_to_absolute(position);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_coordinate_start(peek))
                    break;
            }
        };
        void parse_horizontal_lineto()
        {
            assert(is_horizontal_lineto_start(peek));
            bool is_relative = get() == 'h';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_horizontal_lineto_argument_sequence_start(peek))
                parse_error("missing horizontal lineto argument sequence");
            parse_horizontal_lineto_argument_sequence(is_relative);
        };
        void parse_lineto_argument_sequence(bool is_relative)
        {
            assert(is_lineto_argument_sequence_start(peek));
            while(true)
            {
                auto position = parse_coordinate_pair();
                if(is_relative)
                    path_interpreter.line_to_relative(position);
                else
                    path_interpreter.line_to_absolute(position);
                if(is_comma_wsp_start(peek))
                    parse_comma_wsp();
                if(!is_coordinate_start(peek))
                    break;
            }
        };
        void parse_lineto()
        {
            assert(is_lineto_start(peek));
            bool is_relative = get() == 'l';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_lineto_argument_sequence_start(peek))
                parse_error("missing lineto argument sequence");
            parse_lineto_argument_sequence(is_relative);
        };
        void parse_closepath()
        {
            assert(is_closepath_start(peek));
            next();
            path_interpreter.close_path();
        };
        void parse_moveto_argument_sequence(bool is_relative)
        {
            assert(is_moveto_argument_sequence_start(peek));
            auto position = parse_coordinate_pair();
            if(is_relative)
                path_interpreter.move_to_relative(position);
            else
                path_interpreter.move_to_absolute(position);
            if(is_comma_wsp_start(peek))
                parse_comma_wsp();
            if(is_lineto_argument_sequence_start(peek))
                parse_lineto_argument_sequence(is_relative);
        };
        void parse_moveto()
        {
            assert(is_moveto_start(peek));
            bool is_relative = get() == 'm';
            while(is_wsp(peek))
                parse_wsp();
            if(!is_moveto_argument_sequence_start(peek))
                parse_error("missing moveto argument sequence");
            parse_moveto_argument_sequence(is_relative);
        };
        void parse_drawto_command()
        {
            assert(is_drawto_command_start(peek));
            if(is_closepath_start(peek))
            {
                parse_closepath();
            }
            else if(is_lineto_start(peek))
            {
                parse_lineto();
            }
            else if(is_horizontal_lineto_start(peek))
            {
                parse_horizontal_lineto();
            }
            else if(is_vertical_lineto_start(peek))
            {
                parse_vertical_lineto();
            }
            else if(is_curveto_start(peek))
            {
                parse_curveto();
            }
            else if(is_smooth_curveto_start(peek))
            {
                parse_smooth_curveto();
            }
            else if(is_quadratic_bezier_curveto_start(peek))
            {
                parse_quadratic_bezier_curveto();
            }
            else if(is_smooth_quadratic_bezier_curveto_start(peek))
            {
                parse_smooth_quadratic_bezier_curveto();
            }
            else
            {
                assert(is_elliptical_arc_start(peek));
                parse_elliptical_arc();
            }
        };
        void parse_drawto_commands()
        {
            assert(is_drawto_commands_start(peek));
            while(true)
            {
                parse_drawto_command();
                while(is_wsp(peek))
                    parse_wsp();
                if(!is_drawto_command_start(peek))
                    break;
            }
        };
        void parse_moveto_drawto_command_group()
        {
            assert(is_moveto_drawto_command_group_start(peek));
            parse_moveto();
            while(is_wsp(peek))
                parse_wsp();
            if(is_drawto_commands_start(peek))
                parse_drawto_commands();
        };
        void parse_moveto_drawto_command_groups()
        {
            assert(is_moveto_drawto_command_groups_start(peek));
            while(true)
            {
                parse_moveto_drawto_command_group();
                while(is_wsp(peek))
                    parse_wsp();
                if(!is_moveto_drawto_command_group_start(peek))
                    break;
            }
        };
        void parse()
        {
            next();
            while(is_wsp(peek))
                parse_wsp();
            if(is_moveto_drawto_command_groups_start(peek))
                parse_moveto_drawto_command_groups();
            while(is_wsp(peek))
                parse_wsp();
            if(peek != eof)
                parse_error("unexpected character");
        }
    };
    Parser parser{
        .source = source, .source_size = source_size, .source_name = source_name,
    };
    parser.parse();
    return std::move(parser.path_interpreter).finish();
}

void write_svg_path(std::ostream &os, const Cubic_curve *curves, std::size_t curve_count)
{
    Vec2 last_position{};
    bool has_last_position = false;
    bool is_line_mode = true;
    for(std::size_t i = 0; i < curve_count; i++)
    {
        auto &curve = curves[i];
        bool is_line = curve == Cubic_curve(curve.p0, curve.p3);
        if(!has_last_position || curve.p0 != last_position)
        {
            os << "M" << curve.p0.x << "," << curve.p0.y << "\n";
            is_line_mode = true;
        }
        if(is_line)
        {
            if(!is_line_mode)
                os << "L";
            is_line_mode = true;
            os << curve.p3.x << "," << curve.p3.y << "\n";
        }
        else
        {
            if(is_line_mode)
                os << "C";
            is_line_mode = false;
            os << curve.p1.x << "," << curve.p1.y << " " << curve.p2.x << "," << curve.p2.y << " "
               << curve.p3.x << "," << curve.p3.y << "\n";
        }
        has_last_position = true;
        last_position = curve.p3;
    }
}

void write_svg_path(std::ostream &os, const Line *lines, std::size_t line_count)
{
    Vec2 last_position{};
    bool has_last_position = false;
    for(std::size_t i = 0; i < line_count; i++)
    {
        auto &line = lines[i];
        if(!has_last_position || line.p0 != last_position)
            os << "M" << line.p0.x << "," << line.p0.y << "\n";
        os << line.p1.x << "," << line.p1.y << "\n";
        has_last_position = true;
        last_position = line.p1;
    }
}

template <typename Curve_type>
void write_svg(std::ostream &os,
               const Curve_type *curves,
               std::size_t curve_count,
               const std::string &stroke,
               const std::string &fill)
{
    os << R"(<svg>
    <path d=")";
    write_svg_path(os, curves, curve_count);
    os << R"(" stroke=")" << stroke << R"(" fill=")" << fill << R"("/>
</svg>
)";
}

int main()
{
    try
    {
        std::string path = R"(
M 32,0
A 32,32 0 0 0 0,32
32,32 0 0 0 32,64
32,32 0 0 0 64,32
32,32 0 0 0 32,0
Z
M 16,48 q 16,12 32,0
M 16,16 16,20
M 48,16 48,20
)";
        auto curves = parse_svg_path(path.data(), path.size(), "builtin");
#if 0
        write_svg(std::cout, curves.data(), curves.size(), "blue", "none");
#else
        auto lines = convert_to_lines(curves.data(), curves.size());
        write_svg(std::cout, lines.data(), lines.size(), "blue", "none");
#endif
    }
    catch(Svg_parse_error &e)
    {
        std::cerr << e.what() << std::endl;
    }
}
