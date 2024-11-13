# particles.pyx
import math
import random
from kivy.metrics import dp
from kivy.clock import Clock
cimport cython

cdef class Particle:
    # Type definitions for attributes
    cdef public double base_x, base_y, x, y
    cdef public double size, initial_size
    cdef public tuple color, initial_color
    cdef public double lifetime, age, fade_start_time, fade_duration
    cdef public tuple velocity, acceleration
    cdef public double rotation, rotation_speed
    cdef public double growth_factor, velocity_decay
    cdef public double sin_amplitude_x, sin_frequency_x, sin_phase_x
    cdef public double sin_amplitude_y, sin_frequency_y, sin_phase_y
    cdef public str shape

    def __init__(
        self,
        double x,
        double y,
        double size,
        tuple color,
        double lifetime,
        tuple velocity,
        shape="circle",
        double rotation=0,
        double rotation_speed=0,
        fade_start_time=None,
        fade_duration=None,
        tuple acceleration=(0, 0),
        double brightness_variation=0,
        double growth_factor=1.0,
        double velocity_decay=1.0,
        double sin_amplitude_x=0.0,
        double sin_frequency_x=0.0,
        sin_phase_x=None,
        double sin_amplitude_y=0.0,
        double sin_frequency_y=0.0,
        sin_phase_y=None,
    ):
        self.base_x = x
        self.base_y = y
        self.x = x
        self.y = y
        self.size = dp(size)
        self.initial_size = dp(size)
        self.color = tuple(color)
        self.initial_color = tuple(color)
        self.lifetime = lifetime
        self.velocity = tuple(velocity)
        self.rotation = rotation
        self.rotation_speed = rotation_speed
        self.age = 0
        self.fade_start_time = fade_start_time or lifetime
        self.fade_duration = fade_duration or 0
        self.acceleration = tuple(acceleration)
        self.shape = shape
        self.growth_factor = growth_factor
        self.velocity_decay = velocity_decay

        # Sinusoidal movement parameters
        self.sin_amplitude_x = sin_amplitude_x
        self.sin_frequency_x = sin_frequency_x
        self.sin_phase_x = sin_phase_x if sin_phase_x is not None else random.uniform(0, 2 * math.pi)
        self.sin_amplitude_y = sin_amplitude_y
        self.sin_frequency_y = sin_frequency_y
        self.sin_phase_y = sin_phase_y if sin_phase_y is not None else random.uniform(0, 2 * math.pi)

        # Apply brightness variation
        if brightness_variation > 0:
            brightness_adjustment = random.uniform(-brightness_variation, brightness_variation)
            _color = list(self.color)
            _color[0] = min(max(self.color[0] + brightness_adjustment, 0), 1)
            _color[1] = min(max(self.color[1] + brightness_adjustment, 0), 1)
            _color[2] = min(max(self.color[2] + brightness_adjustment, 0), 1)
            self.color = tuple(_color)

    cpdef void update(self, double dt):
        # Declare variables used within this function
        cdef double sinusoidal_offset_x = 0.0
        cdef double sinusoidal_offset_y = 0.0
        cdef double age_ratio
        cdef double k
        cdef double fade_progress

        # Update age
        self.age += dt

        # Update velocity based on acceleration
        _velocity = list(self.velocity)
        _velocity[0] += self.acceleration[0] * dt
        _velocity[1] += self.acceleration[1] * dt

        # Apply velocity decay (reduce velocity over time)
        _velocity[0] *= self.velocity_decay
        _velocity[1] *= self.velocity_decay

        # Update base position based on velocity
        self.base_x += _velocity[0] * dt
        self.base_y += _velocity[1] * dt

        self.velocity = tuple(_velocity)

        # Update rotation
        self.rotation += self.rotation_speed * dt

        # Compute sinusoidal offsets
        if self.sin_amplitude_x != 0 and self.sin_frequency_x != 0:
            sinusoidal_offset_x = self.sin_amplitude_x * math.sin(2 * math.pi * self.sin_frequency_x * self.age + self.sin_phase_x)

        if self.sin_amplitude_y != 0 and self.sin_frequency_y != 0:
            sinusoidal_offset_y = self.sin_amplitude_y * math.sin(2 * math.pi * self.sin_frequency_y * self.age + self.sin_phase_y)

        # Update final position with sinusoidal movement
        self.x = self.base_x + sinusoidal_offset_x
        self.y = self.base_y + sinusoidal_offset_y

        # Growth or shrinkage over time
        if self.growth_factor != 1.0:
            age_ratio = self.age / self.lifetime
            if self.growth_factor < 1.0:
                k = -math.log(self.growth_factor)
                self.size = self.initial_size * math.exp(-k * age_ratio)
            else:
                k = math.log(self.growth_factor)
                self.size = self.initial_size * math.exp(k * age_ratio)

        # Fading
        if self.age >= self.fade_start_time:
            fade_progress = min(1.0, (self.age - self.fade_start_time) / self.fade_duration)
            _color = list(self.color)
            _color[3] = self.initial_color[3] * (1.0 - fade_progress)
            self.color = tuple(_color)

    cpdef bint is_alive(self):
        return self.age < self.lifetime


cdef class ParticleEmitter:
    # Type definitions for attributes
    cdef public double x, y
    cdef int particle_count
    cdef double particle_lifetime, particle_size
    cdef tuple color
    cdef tuple gravity, acceleration
    cdef double spread
    cdef double fade_start_time, fade_duration
    cdef double size_variation
    cdef tuple initial_velocity
    cdef double velocity_variation
    cdef tuple rotation_speed_range
    cdef str shape
    cdef double brightness_variation
    cdef tuple spawn_area
    cdef double growth_factor, velocity_decay
    cdef double sin_amplitude_x, sin_frequency_x, sin_amplitude_y, sin_frequency_y
    cdef int delay
    cdef list particles
    cdef double time_since_last_emission
    cdef double start_time
    cdef public bint active
    cdef public bint finished
    cdef double emission_rate

    def __init__(
        self,
        double x,
        double y,
        int particle_count,
        double particle_lifetime,
        double particle_size,
        tuple color,
        tuple gravity=(0, -9.8),
        double spread=math.pi / 4,
        double emission_rate=0,
        fade_start_time=0,
        fade_duration=0,
        double size_variation=0.5,
        tuple initial_velocity=(0, 0),
        double velocity_variation=5,
        tuple acceleration=(0, 0),
        tuple rotation_speed_range=(-math.pi, math.pi),
        int delay=0,
        shape="circle",
        double brightness_variation=0,
        tuple spawn_area=(0, 0),
        double growth_factor=1.0,
        double velocity_decay=1.0,
        double sin_amplitude_x=0.0,
        double sin_frequency_x=0.0,
        double sin_amplitude_y=0.0,
        double sin_frequency_y=0.0,
    ):
        self.x = x
        self.y = y
        self.particle_count = particle_count
        self.particle_lifetime = particle_lifetime
        self.particle_size = particle_size
        self.color = color
        self.gravity = gravity
        self.spread = spread
        self.emission_rate = emission_rate
        self.fade_start_time = fade_start_time if fade_start_time is not None else particle_lifetime
        self.fade_duration = fade_duration if fade_duration is not None else 0.0
        self.size_variation = size_variation
        self.initial_velocity = initial_velocity
        self.velocity_variation = velocity_variation
        self.acceleration = acceleration
        self.rotation_speed_range = rotation_speed_range
        self.shape = shape
        self.brightness_variation = brightness_variation
        self.spawn_area = spawn_area
        self.growth_factor = growth_factor
        self.velocity_decay = velocity_decay
        self.sin_amplitude_x = sin_amplitude_x
        self.sin_frequency_x = sin_frequency_x
        self.sin_amplitude_y = sin_amplitude_y
        self.sin_frequency_y = sin_frequency_y
        self.particles = []
        self.active = False
        self.time_since_last_emission = 0.0
        self.delay = delay
        self.start_time = 0
        self.finished = False

    cpdef void emit(self, int count):
        self.active = True
        cdef int particles_to_emit = count if count != 0 else self.particle_count
        for _ in range(particles_to_emit):
            self._create_particle()

    cdef void _create_particle(self):
        cdef double base_angle, angle, base_speed, speed
        cdef tuple velocity
        cdef double lifetime, size, rotation, rotation_speed, spawn_offset_x, spawn_offset_y
        cdef double spawn_x, spawn_y, sin_phase_x, sin_phase_y

        # Calculate the base angle from the initial velocity
        base_angle = math.atan2(self.initial_velocity[1], self.initial_velocity[0])
        angle = base_angle + random.uniform(-self.spread, self.spread)

        # Calculate speed
        base_speed = math.hypot(self.initial_velocity[0], self.initial_velocity[1])
        speed = base_speed + random.uniform(0, self.velocity_variation)

        # Calculate final velocity
        velocity = (speed * math.cos(angle), speed * math.sin(angle))

        # Determine the particle's lifetime and size
        lifetime = random.uniform(0.5 * self.particle_lifetime, self.particle_lifetime)
        size = random.uniform(
            self.particle_size * (1 - self.size_variation),
            self.particle_size * (1 + self.size_variation),
        )
        rotation = random.uniform(0, 2 * math.pi)
        rotation_speed = random.uniform(*self.rotation_speed_range)

        # Determine initial position with spawn area
        spawn_offset_x = random.uniform(-self.spawn_area[0] / 2, self.spawn_area[0] / 2)
        spawn_offset_y = random.uniform(-self.spawn_area[1] / 2, self.spawn_area[1] / 2)
        spawn_x = self.x + spawn_offset_x
        spawn_y = self.y + spawn_offset_y

        # Generate random phase for sinusoidal movement
        sin_phase_x = random.uniform(0, 2 * math.pi)
        sin_phase_y = random.uniform(0, 2 * math.pi)

        particle = Particle(
            spawn_x,
            spawn_y,
            size,
            self.color,
            lifetime,
            velocity,
            shape=self.shape,
            rotation=rotation,
            rotation_speed=rotation_speed,
            fade_start_time=self.fade_start_time,
            fade_duration=self.fade_duration,
            acceleration=(
                self.acceleration[0] + self.gravity[0],
                self.acceleration[1] + self.gravity[1],
            ),
            brightness_variation=self.brightness_variation,
            growth_factor=self.growth_factor,
            velocity_decay=self.velocity_decay,
            sin_amplitude_x=self.sin_amplitude_x,
            sin_frequency_x=self.sin_frequency_x,
            sin_phase_x=sin_phase_x,
            sin_amplitude_y=self.sin_amplitude_y,
            sin_frequency_y=self.sin_frequency_y,
            sin_phase_y=sin_phase_y,
        )
        self.particles.append(particle)

    cpdef void start(self):
        self.start_time = Clock.get_time()
        self.active = True

    cpdef void update(self, double dt):
        if self.start_time is None:
            return

        cdef double current_time = Clock.get_time()
        if current_time - self.start_time < self.delay:
            return

        if not self.active:
            return

        # Continuous emission
        if self.emission_rate:
            self.time_since_last_emission += dt
            while self.time_since_last_emission >= 1 / self.emission_rate:
                self._create_particle()
                self.time_since_last_emission -= 1 / self.emission_rate
        elif not self.particles:  # One-time emission
            self.emit(0)

        # Update particles
        self.particles = [p for p in self.particles if p.is_alive()]
        for particle in self.particles:
            particle.update(dt)

        if not self.particles and not self.emission_rate:
            self.active = False
            self.finished = True

    cpdef list get_particles(self):
        return self.particles
