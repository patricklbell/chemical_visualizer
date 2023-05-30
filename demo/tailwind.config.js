/** @type {import('tailwindcss').Config} */
export default {
  content: ["./index.html", "./src/**/*.{js,ts,jsx,tsx}"],
  darkMode: "class",
  theme: {
    extend: {
      colors: {
        back: "rgb(var(--color-back) / <alpha-value>)",
        back1: "rgb(var(--color-back1) / <alpha-value>)",
        back2: "rgb(var(--color-back2) / <alpha-value>)",
        fore: "rgb(var(--color-fore) / <alpha-value>)",
        fore1: "rgb(var(--color-fore1) / <alpha-value>)",
        fore2: "rgb(var(--color-fore2) / <alpha-value>)",
        accent: "rgb(var(--color-accent) / <alpha-value>)",
        accent1: "rgb(var(--color-accent1) / <alpha-value>)",
        accent2: "rgb(var(--color-accent2) / <alpha-value>)",
        accent3: "rgb(var(--color-accent3) / <alpha-value>)",
        bccent: "rgb(var(--color-bccent) / <alpha-value>)",
        bccent1: "rgb(var(--color-bccent1) / <alpha-value>)",
      },
    },
  },
  plugins: [],
};
