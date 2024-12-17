import { clsx } from "clsx"
import { twMerge } from "tailwind-merge"

export function cn(...inputs) {
  return twMerge(clsx(inputs))
}

export function formatNumber(value) {
  if (typeof value !== 'number') return value;
  
  // For very small numbers, use scientific notation
  if (Math.abs(value) < 0.0001) {
    return value.toExponential(5);
  }
  
  // For regular numbers, use fixed precision
  return Number(value.toPrecision(6)).toString();
}
